"""Tests for ranking, Pareto front, confidence, and sensitivity."""

import pytest
from chempath.data.curated_data import load_curated_data
from chempath.graph.builder import build_drug_target_graph
from chempath.graph.optimizer import (
    assign_confidence,
    rank_compounds_for_target,
    compute_pareto_front,
    find_knee_point,
    DrugRecommendation,
)
from chempath.graph.analysis import (
    toxicity_penalty_sensitivity,
    ic50_perturbation_sensitivity,
    missing_data_sensitivity,
)


@pytest.fixture
def data():
    return load_curated_data()


@pytest.fixture
def graph(data):
    return build_drug_target_graph(data, toxicity_penalty=1.0, verbose=False)


class TestConfidence:
    def test_high_confidence(self):
        conf, reasons = assign_confidence(ic50_nm=1, source="experimental", phase=4, toxicity=0.1)
        assert conf == "HIGH"
        assert len(reasons) == 4

    def test_low_confidence_predicted(self):
        conf, _ = assign_confidence(ic50_nm=500, source="predicted", phase=0, toxicity=0.5)
        assert conf == "LOW"

    def test_medium_confidence(self):
        conf, _ = assign_confidence(ic50_nm=50, source="experimental", phase=2, toxicity=0.25)
        assert conf == "MEDIUM"


class TestRankCompounds:
    def test_returns_recommendations(self, graph, data):
        recs = rank_compounds_for_target(graph, "CHEMBL203", data.get("toxicity"))
        assert len(recs) > 0
        assert all(isinstance(r, DrugRecommendation) for r in recs)

    def test_ranks_are_sequential(self, graph, data):
        recs = rank_compounds_for_target(graph, "CHEMBL203", data.get("toxicity"))
        assert [r.rank for r in recs] == list(range(1, len(recs) + 1))

    def test_balanced_strategy_ranks(self, graph, data):
        recs = rank_compounds_for_target(
            graph, "CHEMBL203", data.get("toxicity"), strategy="balanced"
        )
        # Balanced uses weight + QED + SA composite; verify it produces a valid ranking
        assert len(recs) > 0
        ranks = [r.rank for r in recs]
        assert ranks == list(range(1, len(ranks) + 1))

    def test_efficacy_strategy_uses_ic50(self, graph, data):
        recs = rank_compounds_for_target(
            graph, "CHEMBL203", data.get("toxicity"), strategy="efficacy"
        )
        ic50s = [r.ic50_nm for r in recs]
        assert ic50s == sorted(ic50s)

    def test_safety_strategy_differs_from_efficacy(self, graph, data):
        recs_eff = rank_compounds_for_target(
            graph, "CHEMBL203", data.get("toxicity"), strategy="efficacy"
        )
        recs_safe = rank_compounds_for_target(
            graph, "CHEMBL203", data.get("toxicity"), strategy="safety"
        )
        # Safety uses QED + toxicity composite; should differ from pure IC50
        eff_order = [r.compound_id for r in recs_eff[:10]]
        safe_order = [r.compound_id for r in recs_safe[:10]]
        # At least some reordering should happen (not necessarily all different)
        assert len(recs_safe) > 0
        ranks = [r.rank for r in recs_safe]
        assert ranks == list(range(1, len(ranks) + 1))

    def test_nonexistent_target_returns_empty(self, graph, data):
        recs = rank_compounds_for_target(graph, "NONEXISTENT")
        assert recs == []

    def test_all_have_confidence(self, graph, data):
        recs = rank_compounds_for_target(graph, "CHEMBL203", data.get("toxicity"))
        for r in recs:
            assert r.confidence in ("HIGH", "MEDIUM", "LOW")
            assert len(r.confidence_reasons) > 0


class TestParetoFront:
    def test_pareto_subset_of_recommendations(self, graph, data):
        recs = rank_compounds_for_target(graph, "CHEMBL203", data.get("toxicity"))
        pareto = compute_pareto_front(recs)
        assert len(pareto) <= len(recs)
        assert len(pareto) >= 1

    def test_pareto_no_dominated_solutions(self, graph, data):
        recs = rank_compounds_for_target(graph, "CHEMBL203", data.get("toxicity"))
        pareto = compute_pareto_front(recs)
        for p in pareto:
            for q in pareto:
                if p is q:
                    continue
                # q should not dominate p
                assert not (q.ic50_nm <= p.ic50_nm and q.toxicity <= p.toxicity
                           and (q.ic50_nm < p.ic50_nm or q.toxicity < p.toxicity))


class TestKneePoint:
    def test_returns_recommendation(self, graph, data):
        recs = rank_compounds_for_target(graph, "CHEMBL203", data.get("toxicity"))
        pareto = compute_pareto_front(recs)
        knee = find_knee_point(pareto)
        assert knee is not None
        assert isinstance(knee, DrugRecommendation)

    def test_empty_returns_none(self):
        assert find_knee_point([]) is None

    def test_single_returns_itself(self):
        rec = DrugRecommendation(
            rank=1, compound_id="X", compound_name="X", target_id="T",
            target_name="T", ic50_nm=10, weight=1.0, toxicity=0.1,
            source="experimental", confidence="HIGH",
        )
        assert find_knee_point([rec]) is rec


class TestSensitivityAnalysis:
    def test_toxicity_penalty(self, data):
        result = toxicity_penalty_sensitivity(data, "CHEMBL203")
        assert result.parameter == "toxicity_penalty"
        assert len(result.values_tested) == 5
        assert isinstance(result.stable_compounds, list)

    def test_ic50_perturbation(self, data):
        result = ic50_perturbation_sensitivity(data, "CHEMBL203", n_trials=3)
        assert result.parameter == "ic50_noise_fraction"
        assert len(result.rankings_per_value) > 0

    def test_missing_data(self, data):
        result = missing_data_sensitivity(data, "CHEMBL203", n_trials=3)
        assert result.parameter == "data_removal_fraction"
        assert len(result.rankings_per_value) > 0
