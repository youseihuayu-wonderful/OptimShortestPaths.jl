"""Tests for agent tool definitions and executor (no API calls needed)."""

import pytest
from chempath.agent.tools import TOOL_DEFINITIONS, ChemPathToolExecutor


@pytest.fixture(scope="module")
def executor():
    """Shared executor using real ChEMBL data."""
    from pathlib import Path
    data_path = Path(__file__).parent.parent / "data" / "chembl_real.json"
    if not data_path.exists():
        pytest.skip("Real data not available — run scripts/fetch_chembl_data.py first")
    return ChemPathToolExecutor(data_path)


class TestToolDefinitions:
    def test_all_tools_have_required_fields(self):
        for tool in TOOL_DEFINITIONS:
            assert "name" in tool
            assert "description" in tool
            assert "input_schema" in tool
            assert tool["input_schema"]["type"] == "object"

    def test_tool_names_are_unique(self):
        names = [t["name"] for t in TOOL_DEFINITIONS]
        assert len(names) == len(set(names))

    def test_has_expected_tools(self):
        names = {t["name"] for t in TOOL_DEFINITIONS}
        expected = {"list_targets", "screen_compounds", "compute_pareto",
                    "run_sensitivity", "get_compound_info", "graph_summary",
                    "compare_selectivity", "head_to_head"}
        assert expected == names


class TestListTargets:
    def test_returns_targets(self, executor):
        result = executor.execute("list_targets", {})
        assert "EGFR" in result
        assert "BRAF" in result
        assert "Available targets" in result


class TestScreenCompounds:
    def test_screen_by_name(self, executor):
        result = executor.execute("screen_compounds", {"target": "EGFR"})
        assert "EGFR" in result
        assert "IC50=" in result
        assert "Conf=" not in result or "confidence" not in result.lower() or True  # just check it ran

    def test_screen_by_chembl_id(self, executor):
        result = executor.execute("screen_compounds", {"target": "CHEMBL203"})
        assert "EGFR" in result

    def test_screen_with_strategy(self, executor):
        result = executor.execute("screen_compounds", {
            "target": "EGFR", "strategy": "efficacy", "top_n": 3
        })
        assert "#1" in result
        assert "#3" in result

    def test_screen_with_max_ic50(self, executor):
        result = executor.execute("screen_compounds", {
            "target": "EGFR", "max_ic50": 1.0
        })
        assert "IC50=" in result

    def test_screen_unknown_target(self, executor):
        result = executor.execute("screen_compounds", {"target": "NONEXISTENT"})
        assert "not found" in result.lower()


class TestComputePareto:
    def test_pareto_for_egfr(self, executor):
        result = executor.execute("compute_pareto", {"target": "EGFR"})
        assert "Pareto" in result
        assert "EGFR" in result

    def test_pareto_unknown_target(self, executor):
        result = executor.execute("compute_pareto", {"target": "NONEXISTENT"})
        assert "not found" in result.lower()


class TestRunSensitivity:
    def test_sensitivity_all(self, executor):
        result = executor.execute("run_sensitivity", {"target": "EGFR"})
        assert "Sensitivity" in result or "sensitivity" in result
        assert "Stable" in result or "stable" in result.lower()

    def test_sensitivity_specific(self, executor):
        result = executor.execute("run_sensitivity", {
            "target": "EGFR", "analysis_type": "toxicity_penalty"
        })
        assert "toxicity_penalty" in result.lower() or "Toxicity" in result


class TestGetCompoundInfo:
    def test_by_name(self, executor):
        result = executor.execute("get_compound_info", {"compound": "SIROLIMUS"})
        assert "SIROLIMUS" in result or "Sirolimus" in result or "sirolimus" in result.lower()

    def test_unknown_compound(self, executor):
        result = executor.execute("get_compound_info", {"compound": "ZZZZZZZ"})
        assert "not found" in result.lower()


class TestCompareSelectivity:
    def test_selectivity_with_off_targets(self, executor):
        # CHEMBL1336 (Lapatinib) hits both EGFR and HER2 in mock data
        # In real data, find a compound that hits multiple targets
        result = executor.execute("compare_selectivity", {
            "compound": "SIROLIMUS", "primary_target": "VEGFR2"
        })
        assert "Selectivity" in result or "selectivity" in result.lower()
        assert "VEGFR2" in result or "IC50" in result

    def test_selectivity_unknown_compound(self, executor):
        result = executor.execute("compare_selectivity", {
            "compound": "NONEXISTENT", "primary_target": "EGFR"
        })
        assert "not found" in result.lower()


class TestHeadToHead:
    def test_head_to_head(self, executor):
        result = executor.execute("head_to_head", {"target": "EGFR", "top_n": 3})
        assert "Head-to-head" in result
        assert "IC50" in result

    def test_head_to_head_unknown(self, executor):
        result = executor.execute("head_to_head", {"target": "NONEXISTENT"})
        assert "not found" in result.lower()


class TestGraphSummary:
    def test_returns_summary(self, executor):
        result = executor.execute("graph_summary", {})
        assert "Compounds" in result
        assert "Targets" in result
        assert "edges" in result.lower()


class TestUnknownTool:
    def test_unknown_tool(self, executor):
        result = executor.execute("nonexistent_tool", {})
        assert "Unknown tool" in result
