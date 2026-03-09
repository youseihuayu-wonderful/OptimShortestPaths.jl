"""Tests for graph building and summary."""

import math
import pytest
from chempath.data.curated_data import load_curated_data
from chempath.graph.builder import (
    ic50_to_weight,
    build_drug_target_graph,
    add_predicted_edges,
    get_graph_summary,
)


class TestIC50ToWeight:
    def test_high_potency(self):
        # 1 nM → pIC50 = 9, weight = 14 - 9 = 5.0
        w = ic50_to_weight(1.0)
        assert w == pytest.approx(5.0, abs=0.01)

    def test_moderate_potency(self):
        # 100 nM → pIC50 = 7, weight = 14 - 7 = 7.0
        w = ic50_to_weight(100.0)
        assert w == pytest.approx(7.0, abs=0.01)

    def test_low_potency(self):
        # 10000 nM → pIC50 = 5, weight = 14 - 5 = 9.0
        w = ic50_to_weight(10000.0)
        assert w == pytest.approx(9.0, abs=0.01)

    def test_zero_returns_inf(self):
        assert ic50_to_weight(0) == float("inf")

    def test_negative_returns_inf(self):
        assert ic50_to_weight(-1) == float("inf")

    def test_monotonic(self):
        """Lower IC50 (more potent) → lower weight."""
        assert ic50_to_weight(1) < ic50_to_weight(10) < ic50_to_weight(100)


class TestBuildGraph:
    @pytest.fixture
    def data(self):
        return load_curated_data()

    @pytest.fixture
    def graph(self, data):
        return build_drug_target_graph(data, toxicity_penalty=1.0, verbose=False)

    def test_node_counts(self, graph):
        summary = get_graph_summary(graph)
        assert summary["compounds"] == 18  # 20 - 2 invalid SMILES
        assert summary["targets"] == 11

    def test_edge_count(self, graph):
        summary = get_graph_summary(graph)
        assert summary["experimental_edges"] == 27
        assert summary["predicted_edges"] == 0

    def test_rejects_invalid_smiles(self, graph):
        assert "CHEMBL_INVALID" not in graph
        assert "CHEMBL_EMPTY" not in graph

    def test_edge_has_weight(self, graph):
        edge = graph.edges["CHEMBL1201585", "CHEMBL203"]
        assert "weight" in edge
        assert edge["weight"] > 0
        assert edge["ic50_nm"] == 0.5

    def test_toxicity_penalty_applied(self, data):
        g0 = build_drug_target_graph(data, toxicity_penalty=0.0, verbose=False)
        g5 = build_drug_target_graph(data, toxicity_penalty=5.0, verbose=False)
        # Sorafenib (CHEMBL1336) → BRAF (CHEMBL5145)
        w0 = g0.edges["CHEMBL1336", "CHEMBL5145"]["weight"]
        w5 = g5.edges["CHEMBL1336", "CHEMBL5145"]["weight"]
        assert w5 > w0  # higher penalty → higher weight

    def test_compound_node_attributes(self, graph):
        node = graph.nodes["CHEMBL1201585"]
        assert node["name"] == "Osimertinib"
        assert node["node_type"] == "compound"
        assert node["phase"] == 4

    def test_target_node_attributes(self, graph):
        node = graph.nodes["CHEMBL203"]
        assert node["name"] == "EGFR"
        assert node["node_type"] == "target"


class TestAddPredictedEdges:
    def test_adds_predicted_edge(self):
        data = load_curated_data()
        G = build_drug_target_graph(data, verbose=False)
        predictions = [
            {"compound": "CHEMBL1201585", "target": "CHEMBL4282", "predicted_ic50": 50.0}
        ]
        G = add_predicted_edges(G, predictions, uncertainty_penalty=0.2)
        assert G.has_edge("CHEMBL1201585", "CHEMBL4282")
        edge = G.edges["CHEMBL1201585", "CHEMBL4282"]
        assert edge["source"] == "predicted"

    def test_does_not_overwrite_experimental(self):
        data = load_curated_data()
        G = build_drug_target_graph(data, verbose=False)
        original_weight = G.edges["CHEMBL1201585", "CHEMBL203"]["weight"]
        predictions = [
            {"compound": "CHEMBL1201585", "target": "CHEMBL203", "predicted_ic50": 999.0}
        ]
        G = add_predicted_edges(G, predictions)
        assert G.edges["CHEMBL1201585", "CHEMBL203"]["weight"] == original_weight
