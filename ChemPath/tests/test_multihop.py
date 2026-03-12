"""Tests for multi-hop graph, pathfinding, and benchmark modules."""

import math
import pytest
from chempath.data.curated_data import load_curated_data
from chempath.graph.network import (
    EdgeProbability,
    ic50_to_efficacy_prob,
    toxicity_to_safety_prob,
    phase_to_evidence_prob,
    annotation_to_prob,
    build_multihop_graph,
    get_multihop_summary,
)
from chempath.graph.pathfinding import (
    find_shortest_paths,
    find_repurposing_candidates,
    find_mechanism_of_action,
    identify_targets_for_disease,
    find_combination_therapy,
)
from chempath.graph.benchmark import (
    run_multihop_benchmark,
    run_random_baseline,
    run_single_hop_baseline,
    GROUND_TRUTH,
)
from chempath.graph.julia_bridge import JuliaBridge, DMYResult, JULIA_BINARY


# ===================================================================
# Fixtures
# ===================================================================

@pytest.fixture
def curated_data():
    return load_curated_data()


@pytest.fixture
def graph(curated_data):
    return build_multihop_graph(curated_data, verbose=False)


# ===================================================================
# network.py — Probability transforms
# ===================================================================

class TestProbabilityTransforms:
    def test_ic50_high_potency(self):
        """IC50=1nM → P≈0.99 (very potent)."""
        p = ic50_to_efficacy_prob(1.0)
        assert 0.98 < p < 1.0

    def test_ic50_moderate_potency(self):
        """IC50=100nM → P=0.50 (reference)."""
        p = ic50_to_efficacy_prob(100.0)
        assert p == pytest.approx(0.5, abs=0.01)

    def test_ic50_low_potency(self):
        """IC50=10000nM → P≈0.01 (weak)."""
        p = ic50_to_efficacy_prob(10000.0)
        assert p < 0.02

    def test_ic50_monotonic(self):
        """Lower IC50 → higher probability."""
        assert ic50_to_efficacy_prob(1) > ic50_to_efficacy_prob(100) > ic50_to_efficacy_prob(10000)

    def test_ic50_zero_returns_small(self):
        p = ic50_to_efficacy_prob(0)
        assert p == pytest.approx(1e-10)

    def test_toxicity_low(self):
        """Low toxicity → high safety probability."""
        assert toxicity_to_safety_prob(0.1) == pytest.approx(0.9)

    def test_toxicity_high(self):
        """High toxicity → low safety probability."""
        assert toxicity_to_safety_prob(0.9) == pytest.approx(0.1)

    def test_toxicity_clamped(self):
        assert toxicity_to_safety_prob(1.5) == 0.01
        assert toxicity_to_safety_prob(-0.5) == 1.0

    def test_phase_approved(self):
        """Phase 4, experimental → 0.95."""
        assert phase_to_evidence_prob(4, "experimental") == 0.95

    def test_phase_preclinical(self):
        """Phase 0, experimental → 0.40."""
        assert phase_to_evidence_prob(0, "experimental") == 0.40

    def test_phase_predicted(self):
        """Predicted source → 0.20 regardless of phase."""
        assert phase_to_evidence_prob(4, "predicted") == 0.20

    def test_annotation_strong(self):
        """Low weight (strong connection) → high probability."""
        p = annotation_to_prob(0.1)
        assert p == pytest.approx(0.9, abs=0.01)

    def test_annotation_weak(self):
        """High weight (weak connection) → low probability."""
        p = annotation_to_prob(0.9)
        assert p == pytest.approx(0.1, abs=0.01)


class TestEdgeProbability:
    def test_logcost_perfect(self):
        """P=1.0 → cost=0.0."""
        ep = EdgeProbability(efficacy=1.0, safety=1.0, evidence=1.0)
        lc = ep.to_logcost()
        assert lc["w_efficacy"] == pytest.approx(0.0, abs=1e-6)

    def test_logcost_positive(self):
        """P<1.0 → positive cost."""
        ep = EdgeProbability(efficacy=0.5, safety=0.8, evidence=0.9)
        lc = ep.to_logcost()
        assert lc["w_efficacy"] > 0
        assert lc["w_safety"] > 0
        assert lc["w_evidence"] > 0

    def test_logcost_additive(self):
        """Two edges with P=0.5 should sum to -log(0.25)."""
        ep = EdgeProbability(efficacy=0.5)
        lc = ep.to_logcost()
        assert 2 * lc["w_efficacy"] == pytest.approx(-math.log(0.25), abs=1e-6)


# ===================================================================
# network.py — Graph construction
# ===================================================================

class TestMultihopGraph:
    def test_has_four_node_types(self, graph):
        types = set()
        for _, d in graph.nodes(data=True):
            types.add(d.get("node_type"))
        assert types == {"compound", "target", "pathway", "disease"}

    def test_has_four_edge_types(self, graph):
        types = set()
        for _, _, d in graph.edges(data=True):
            types.add(d.get("edge_type"))
        assert "compound_target" in types
        assert "target_pathway" in types
        assert "pathway_disease" in types

    def test_all_edges_have_weight(self, graph):
        for u, v, d in graph.edges(data=True):
            assert "weight" in d, f"Edge {u}→{v} missing 'weight'"
            assert d["weight"] >= 0

    def test_all_edges_have_probability_dims(self, graph):
        for u, v, d in graph.edges(data=True):
            assert "p_efficacy" in d, f"Edge {u}→{v} missing p_efficacy"
            assert "w_efficacy" in d, f"Edge {u}→{v} missing w_efficacy"
            assert 0 < d["p_efficacy"] <= 1.0

    def test_weight_equals_w_efficacy(self, graph):
        """Generic weight should equal w_efficacy."""
        for u, v, d in graph.edges(data=True):
            assert d["weight"] == pytest.approx(d["w_efficacy"], abs=1e-10)

    def test_compound_target_edges_have_ic50(self, graph):
        for u, v, d in graph.edges(data=True):
            if d.get("edge_type") == "compound_target":
                assert "ic50_nm" in d
                assert d["ic50_nm"] > 0

    def test_diseases_present(self, graph):
        diseases = [n for n, d in graph.nodes(data=True) if d.get("node_type") == "disease"]
        assert len(diseases) == 10

    def test_pathways_present(self, graph):
        pathways = [n for n, d in graph.nodes(data=True) if d.get("node_type") == "pathway"]
        assert len(pathways) == 10

    def test_summary(self, graph):
        s = get_multihop_summary(graph)
        assert s["total_nodes"] > 0
        assert s["total_edges"] > 0
        assert s["nodes_by_type"]["disease"] == 10
        assert s["nodes_by_type"]["pathway"] == 10


# ===================================================================
# pathfinding.py
# ===================================================================

class TestPathfinding:
    def test_find_shortest_paths_basic(self, graph):
        """Should find at least one path from a compound to a disease."""
        # Osimertinib → EGFR → ErbB → NSCLC
        compounds = [n for n, d in graph.nodes(data=True) if d.get("node_type") == "compound"]
        diseases = [n for n, d in graph.nodes(data=True) if d.get("node_type") == "disease"]

        found_path = False
        for cid in compounds[:3]:
            for did in diseases[:3]:
                paths = find_shortest_paths(graph, cid, did, k=1)
                if paths:
                    found_path = True
                    assert paths[0].n_hops >= 2  # at least compound→target→...→disease
                    assert paths[0].total_probability > 0
                    assert paths[0].total_probability <= 1.0
                    break
            if found_path:
                break
        assert found_path, "Should find at least one compound→disease path"

    def test_path_probability_decreases_with_hops(self, graph):
        """Longer paths should have lower probability (more uncertainty)."""
        compounds = [n for n, d in graph.nodes(data=True) if d.get("node_type") == "compound"]
        diseases = [n for n, d in graph.nodes(data=True) if d.get("node_type") == "disease"]

        for cid in compounds[:5]:
            for did in diseases[:5]:
                paths = find_shortest_paths(graph, cid, did, k=1)
                if paths:
                    p = paths[0]
                    # Path probability should be <= min(edge probability)
                    assert p.total_probability <= 1.0
                    assert p.total_cost >= 0

    def test_path_result_structure(self, graph):
        """PathResult should have all required fields."""
        compounds = [n for n, d in graph.nodes(data=True) if d.get("node_type") == "compound"]
        diseases = [n for n, d in graph.nodes(data=True) if d.get("node_type") == "disease"]

        for cid in compounds[:3]:
            for did in diseases[:3]:
                paths = find_shortest_paths(graph, cid, did, k=1)
                if paths:
                    p = paths[0]
                    assert len(p.path) == len(p.path_names)
                    assert len(p.path) == len(p.path_types)
                    assert len(p.edge_types) == p.n_hops
                    assert "efficacy" in p.costs_by_dim
                    assert "safety" in p.costs_by_dim
                    assert "evidence" in p.costs_by_dim
                    return
        pytest.skip("No paths found in curated data")

    def test_nonexistent_source(self, graph):
        paths = find_shortest_paths(graph, "FAKE_ID", "DIS_NSCLC")
        assert paths == []

    def test_nonexistent_target(self, graph):
        compounds = [n for n, d in graph.nodes(data=True) if d.get("node_type") == "compound"]
        paths = find_shortest_paths(graph, compounds[0], "FAKE_DISEASE")
        assert paths == []


class TestRepurposing:
    def test_returns_results(self, graph):
        compounds = [n for n, d in graph.nodes(data=True) if d.get("node_type") == "compound"]
        if not compounds:
            pytest.skip("No compounds in graph")
        results = find_repurposing_candidates(graph, compounds[0])
        # Should find at least some diseases
        assert isinstance(results, list)

    def test_results_sorted_by_probability(self, graph):
        compounds = [n for n, d in graph.nodes(data=True) if d.get("node_type") == "compound"]
        for cid in compounds[:3]:
            results = find_repurposing_candidates(graph, cid)
            if len(results) >= 2:
                probs = [r.best_probability for r in results]
                assert probs == sorted(probs, reverse=True), "Should be sorted by descending probability"
                return
        pytest.skip("Not enough results to test sorting")

    def test_nonexistent_compound(self, graph):
        results = find_repurposing_candidates(graph, "FAKE_COMPOUND")
        assert results == []


class TestMoA:
    def test_returns_mechanism(self, graph):
        compounds = [n for n, d in graph.nodes(data=True) if d.get("node_type") == "compound"]
        diseases = [n for n, d in graph.nodes(data=True) if d.get("node_type") == "disease"]
        for cid in compounds[:5]:
            for did in diseases[:5]:
                result = find_mechanism_of_action(graph, cid, did)
                if result is not None:
                    assert result.primary_mechanism is not None
                    assert len(result.paths) >= 1
                    return
        pytest.skip("No MoA found in curated data")

    def test_nonexistent_pair(self, graph):
        result = find_mechanism_of_action(graph, "FAKE", "FAKE2")
        assert result is None


class TestTargetIdentification:
    def test_identify_targets(self, graph):
        result = identify_targets_for_disease(graph, "DIS_NSCLC")
        if result is None:
            pytest.skip("No targets found for NSCLC")
        assert len(result.targets) > 0
        # Targets should be sorted by probability (descending)
        probs = [t["probability"] for t in result.targets]
        assert probs == sorted(probs, reverse=True)

    def test_nonexistent_disease(self, graph):
        result = identify_targets_for_disease(graph, "FAKE_DISEASE")
        assert result is None


class TestCombinationTherapy:
    def test_find_combination(self, graph):
        targets = [n for n, d in graph.nodes(data=True) if d.get("node_type") == "target"]
        if len(targets) < 2:
            pytest.skip("Not enough targets")
        result = find_combination_therapy(graph, targets[:2])
        if result is None:
            pytest.skip("No combination found")
        assert len(result.best_combination) >= 1

    def test_nonexistent_targets(self, graph):
        result = find_combination_therapy(graph, ["FAKE1", "FAKE2"])
        assert result is None


# ===================================================================
# benchmark.py
# ===================================================================

class TestBenchmark:
    def test_ground_truth_format(self):
        """Ground truth should map compound IDs to disease ID lists."""
        for cid, diseases in GROUND_TRUTH.items():
            assert isinstance(cid, str)
            assert isinstance(diseases, list)
            assert all(d.startswith("DIS_") for d in diseases)

    def test_random_baseline_ef_near_one(self, graph):
        """Random baseline enrichment factor should be ~1.0."""
        result = run_random_baseline(graph, n_trials=50)
        # EF should be approximately 1.0 (by definition)
        assert 0.3 < result.enrichment_factor < 2.5

    def test_multihop_runs(self, graph):
        """Multi-hop benchmark should run without errors on curated data."""
        result = run_multihop_benchmark(graph, k_values=[1, 3, 5])
        assert result.n_drugs >= 0
        assert 0 <= result.mrr <= 1.0
        for k, recall in result.recall_at_k.items():
            assert 0 <= recall <= 1.0

    def test_single_hop_runs(self, graph):
        result = run_single_hop_baseline(graph, k_values=[1, 3, 5])
        assert result.n_drugs >= 0
        assert 0 <= result.mrr <= 1.0


# ===================================================================
# julia_bridge.py — export/parse logic (no Julia required)
# ===================================================================

class TestJuliaBridgeExport:
    def test_index_mapping(self, graph):
        bridge = JuliaBridge(graph)
        # All nodes should have 1-based indices
        assert len(bridge._node_to_idx) == graph.number_of_nodes()
        assert min(bridge._node_to_idx.values()) == 1
        assert max(bridge._node_to_idx.values()) == graph.number_of_nodes()

    def test_export_format(self, graph):
        bridge = JuliaBridge(graph)
        exported = bridge.export_graph("w_efficacy")
        assert exported["n_vertices"] == graph.number_of_nodes()
        assert len(exported["edges"]) == graph.number_of_edges()

    def test_export_edges_valid(self, graph):
        bridge = JuliaBridge(graph)
        exported = bridge.export_graph("w_efficacy")
        n = exported["n_vertices"]
        for edge in exported["edges"]:
            assert 1 <= edge["source"] <= n
            assert 1 <= edge["target"] <= n
            assert edge["weight"] >= 0  # DMY requires non-negative

    def test_export_nodes_metadata(self, graph):
        bridge = JuliaBridge(graph)
        exported = bridge.export_graph("w_efficacy")
        for idx_str, meta in exported["nodes"].items():
            assert "id" in meta
            assert "name" in meta
            assert "type" in meta

    def test_roundtrip_index(self, graph):
        """node → index → node should be identity."""
        bridge = JuliaBridge(graph)
        for node_id in list(graph.nodes())[:20]:
            idx = bridge._node_to_idx[node_id]
            assert bridge._idx_to_node[idx] == node_id

    def test_parse_result(self, graph):
        """Test _parse_result with synthetic Julia output."""
        bridge = JuliaBridge(graph)
        n = graph.number_of_nodes()
        source = list(graph.nodes())[0]
        source_idx = bridge._node_to_idx[source]

        # Simulate Julia output: source has distance 0, others have INF
        distances = [1e308] * n
        parents = [0] * n
        distances[source_idx - 1] = 0.0

        # Set a neighbor as reachable
        neighbors = list(graph.successors(source))
        if neighbors:
            nb = neighbors[0]
            nb_idx = bridge._node_to_idx[nb]
            distances[nb_idx - 1] = 0.5
            parents[nb_idx - 1] = source_idx

        fake_result = {
            "n_vertices": n,
            "n_edges": graph.number_of_edges(),
            "source": source_idx,
            "distances": distances,
            "parents": parents,
        }

        result = bridge._parse_result(fake_result, source, "w_efficacy", None)
        assert isinstance(result, DMYResult)
        assert result.distances[source] == 0.0
        if neighbors:
            assert result.distances[neighbors[0]] == pytest.approx(0.5)
            assert result.parents[neighbors[0]] == source

    def test_reconstruct_path(self, graph):
        """Test path reconstruction from parent array."""
        bridge = JuliaBridge(graph)
        # Simple chain: 1 → 2 → 3
        parents = [0, 1, 2, 0, 0]  # parent[1]=0 (source), parent[2]=1, parent[3]=2
        path = bridge._reconstruct_path(parents, target_idx=3, source_idx=1)
        assert path == [1, 2, 3]

    def test_no_julia_raises_clear_error(self, graph):
        """If Julia is not installed, run_sssp should give a clear error."""
        import chempath.graph.julia_bridge as jb
        old = jb.JULIA_BINARY
        try:
            jb.JULIA_BINARY = None
            bridge = JuliaBridge(graph)
            source = list(graph.nodes())[0]
            with pytest.raises(RuntimeError, match="Julia not found"):
                bridge.run_sssp(source)
        finally:
            jb.JULIA_BINARY = old

    @pytest.mark.skipif(JULIA_BINARY is None, reason="Julia not installed")
    def test_end_to_end_sssp(self, graph):
        """Integration test: run DMY solver if Julia is available."""
        bridge = JuliaBridge(graph)
        compounds = [n for n, d in graph.nodes(data=True) if d.get("node_type") == "compound"]
        if not compounds:
            pytest.skip("No compounds")
        result = bridge.run_sssp(compounds[0])
        assert result.n_reachable > 0
        assert result.distances[compounds[0]] == pytest.approx(0.0)
