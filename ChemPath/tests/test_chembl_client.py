"""Tests for ChEMBL client — unit tests with no network calls."""

import pytest
from chempath.data.chembl_client import clean_activity, ANTICANCER_TARGETS


class TestCleanActivity:
    def test_valid_record_nm(self):
        raw = {
            "molecule_chembl_id": "CHEMBL1201585",
            "molecule_pref_name": "OSIMERTINIB",
            "canonical_smiles": "C=CC(=O)Nc1ccccc1",
            "standard_value": "0.5",
            "standard_units": "nM",
            "standard_type": "IC50",
            "target_chembl_id": "CHEMBL203",
            "target_pref_name": "EGFR",
            "pchembl_value": "9.30",
            "assay_type": "B",
            "data_validity_comment": None,
        }
        result = clean_activity(raw)
        assert result is not None
        assert result["compound"] == "CHEMBL1201585"
        assert result["value"] == 0.5
        assert result["units"] == "nM"
        assert result["source"] == "experimental"

    def test_valid_record_um(self):
        raw = {
            "molecule_chembl_id": "CHEMBL100",
            "canonical_smiles": "CCO",
            "standard_value": "1.5",
            "standard_units": "uM",
            "target_chembl_id": "CHEMBL203",
            "data_validity_comment": None,
        }
        result = clean_activity(raw)
        assert result is not None
        assert result["value"] == 1500.0  # 1.5 uM → 1500 nM

    def test_valid_record_pm(self):
        raw = {
            "molecule_chembl_id": "CHEMBL100",
            "canonical_smiles": "CCO",
            "standard_value": "500",
            "standard_units": "pM",
            "target_chembl_id": "CHEMBL203",
            "data_validity_comment": None,
        }
        result = clean_activity(raw)
        assert result is not None
        assert result["value"] == 0.5  # 500 pM → 0.5 nM

    def test_missing_value_returns_none(self):
        raw = {
            "molecule_chembl_id": "CHEMBL100",
            "canonical_smiles": "CCO",
            "standard_value": None,
            "standard_units": "nM",
            "target_chembl_id": "CHEMBL203",
        }
        assert clean_activity(raw) is None

    def test_missing_smiles_returns_none(self):
        raw = {
            "molecule_chembl_id": "CHEMBL100",
            "canonical_smiles": None,
            "standard_value": "10",
            "standard_units": "nM",
            "target_chembl_id": "CHEMBL203",
        }
        assert clean_activity(raw) is None

    def test_zero_value_returns_none(self):
        raw = {
            "molecule_chembl_id": "CHEMBL100",
            "canonical_smiles": "CCO",
            "standard_value": "0",
            "standard_units": "nM",
            "target_chembl_id": "CHEMBL203",
            "data_validity_comment": None,
        }
        assert clean_activity(raw) is None

    def test_negative_value_returns_none(self):
        raw = {
            "molecule_chembl_id": "CHEMBL100",
            "canonical_smiles": "CCO",
            "standard_value": "-5",
            "standard_units": "nM",
            "target_chembl_id": "CHEMBL203",
            "data_validity_comment": None,
        }
        assert clean_activity(raw) is None

    def test_flagged_data_returns_none(self):
        raw = {
            "molecule_chembl_id": "CHEMBL100",
            "canonical_smiles": "CCO",
            "standard_value": "10",
            "standard_units": "nM",
            "target_chembl_id": "CHEMBL203",
            "data_validity_comment": "Outside typical range",
        }
        assert clean_activity(raw) is None

    def test_unknown_units_returns_none(self):
        raw = {
            "molecule_chembl_id": "CHEMBL100",
            "canonical_smiles": "CCO",
            "standard_value": "10",
            "standard_units": "mg/mL",
            "target_chembl_id": "CHEMBL203",
            "data_validity_comment": None,
        }
        assert clean_activity(raw) is None

    def test_non_numeric_value_returns_none(self):
        raw = {
            "molecule_chembl_id": "CHEMBL100",
            "canonical_smiles": "CCO",
            "standard_value": "active",
            "standard_units": "nM",
            "target_chembl_id": "CHEMBL203",
            "data_validity_comment": None,
        }
        assert clean_activity(raw) is None

    def test_missing_name_uses_id(self):
        raw = {
            "molecule_chembl_id": "CHEMBL100",
            "molecule_pref_name": None,
            "canonical_smiles": "CCO",
            "standard_value": "10",
            "standard_units": "nM",
            "target_chembl_id": "CHEMBL203",
            "target_pref_name": None,
            "data_validity_comment": None,
        }
        result = clean_activity(raw)
        assert result["compound_name"] == "CHEMBL100"
        assert result["target_name"] == "CHEMBL203"


class TestTargetList:
    def test_has_ten_targets(self):
        assert len(ANTICANCER_TARGETS) == 10

    def test_all_have_required_fields(self):
        for t in ANTICANCER_TARGETS:
            assert "chembl_id" in t
            assert "name" in t
            assert t["chembl_id"].startswith("CHEMBL")
