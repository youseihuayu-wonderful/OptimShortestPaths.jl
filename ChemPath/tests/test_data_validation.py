"""Tests for automated data validation."""

import copy
import pytest
from chempath.data.curated_data import load_curated_data
from chempath.data.data_validation import validate_dataset


class TestCuratedDataPasses:
    def test_no_errors(self):
        errors = validate_dataset()
        assert errors == [], f"Validation errors:\n" + "\n".join(errors)


class TestDuplicateCompound:
    def test_catches_duplicate_compound_id(self):
        data = copy.deepcopy(load_curated_data())
        data["compounds"].append(data["compounds"][0])
        errors = validate_dataset(data)
        assert any("Duplicate compound ID" in e for e in errors)


class TestDuplicateTarget:
    def test_catches_duplicate_target_id(self):
        data = copy.deepcopy(load_curated_data())
        data["targets"].append(data["targets"][0])
        errors = validate_dataset(data)
        assert any("Duplicate target ID" in e for e in errors)


class TestMislabelDetection:
    def test_catches_wrong_compound_name(self):
        data = copy.deepcopy(load_curated_data())
        for c in data["compounds"]:
            if c["chembl_id"] == "CHEMBL1336":
                c["name"] = "Lapatinib"  # wrong — should be Sorafenib
        errors = validate_dataset(data)
        assert any("mislabel" in e.lower() and "CHEMBL1336" in e for e in errors)

    def test_catches_wrong_target_name(self):
        data = copy.deepcopy(load_curated_data())
        for t in data["targets"]:
            if t["chembl_id"] == "CHEMBL203":
                t["name"] = "WRONG"
        errors = validate_dataset(data)
        assert any("mislabel" in e.lower() and "CHEMBL203" in e for e in errors)


class TestInvalidReferences:
    def test_catches_unknown_compound_in_bioactivity(self):
        data = copy.deepcopy(load_curated_data())
        data["bioactivities"].append({
            "compound": "CHEMBL_FAKE", "target": "CHEMBL203",
            "type": "IC50", "value": 10.0, "units": "nM", "source": "experimental",
        })
        errors = validate_dataset(data)
        assert any("unknown compound" in e.lower() for e in errors)

    def test_catches_unknown_target_in_bioactivity(self):
        data = copy.deepcopy(load_curated_data())
        data["bioactivities"].append({
            "compound": "CHEMBL939", "target": "CHEMBL_FAKE",
            "type": "IC50", "value": 10.0, "units": "nM", "source": "experimental",
        })
        errors = validate_dataset(data)
        assert any("unknown target" in e.lower() for e in errors)


class TestIC50Range:
    def test_catches_negative_ic50(self):
        data = copy.deepcopy(load_curated_data())
        data["bioactivities"][0]["value"] = -5.0
        errors = validate_dataset(data)
        assert any("non-positive" in e for e in errors)

    def test_catches_extremely_high_ic50(self):
        data = copy.deepcopy(load_curated_data())
        data["bioactivities"][0]["value"] = 999_999.0
        errors = validate_dataset(data)
        assert any("above" in e for e in errors)


class TestDuplicateEdges:
    def test_catches_duplicate_bioactivity(self):
        data = copy.deepcopy(load_curated_data())
        data["bioactivities"].append(data["bioactivities"][0])
        errors = validate_dataset(data)
        assert any("duplicate edge" in e.lower() for e in errors)


class TestToxicityCoverage:
    def test_catches_missing_toxicity(self):
        data = copy.deepcopy(load_curated_data())
        del data["toxicity"]["CHEMBL939"]
        errors = validate_dataset(data)
        assert any("CHEMBL939" in e and "no toxicity" in e for e in errors)

    def test_catches_extra_toxicity(self):
        data = copy.deepcopy(load_curated_data())
        data["toxicity"]["CHEMBL_GHOST"] = {"overall": 0.5}
        errors = validate_dataset(data)
        assert any("CHEMBL_GHOST" in e and "unknown compound" in e.lower() for e in errors)
