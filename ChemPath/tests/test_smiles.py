"""Tests for SMILES validation."""

import pytest
from chempath.chemistry.smiles import (
    validate_smiles_basic,
    validate_smiles,
    validate_batch,
    filter_valid_compounds,
    ValidationResult,
)


class TestBasicValidation:
    def test_valid_simple_smiles(self):
        result = validate_smiles_basic("CCO")
        assert result.is_valid is True
        assert result.level == "basic"

    def test_valid_complex_smiles(self):
        # Osimertinib
        result = validate_smiles_basic(
            "C=CC(=O)Nc1cc(OC)c(Nc2nccc(-c3cn(C)c4ccccc34)n2)cc1N(C)CCN(C)C"
        )
        assert result.is_valid is True

    def test_empty_string(self):
        result = validate_smiles_basic("")
        assert result.is_valid is False
        assert "Empty" in result.error

    def test_none_input(self):
        result = validate_smiles_basic(None)
        assert result.is_valid is False

    def test_invalid_characters(self):
        result = validate_smiles_basic("NOT_A_SMILES")
        assert result.is_valid is False
        assert "Invalid characters" in result.error

    def test_unmatched_bracket(self):
        result = validate_smiles_basic("[CH3")
        assert result.is_valid is False
        assert "Unmatched" in result.error

    def test_unmatched_parenthesis(self):
        result = validate_smiles_basic("C(CC")
        assert result.is_valid is False
        assert "Unmatched" in result.error

    def test_no_atoms(self):
        result = validate_smiles_basic("123")
        assert result.is_valid is False
        assert "No atoms" in result.error

    def test_whitespace_stripped(self):
        result = validate_smiles_basic("  CCO  ")
        assert result.is_valid is True
        assert result.smiles == "CCO"


class TestBatchValidation:
    def test_mixed_batch(self):
        smiles = ["CCO", "NOT_VALID", "c1ccccc1"]
        result = validate_batch(smiles)
        assert result["summary"]["total"] == 3
        assert result["summary"]["valid_count"] == 2
        assert result["summary"]["invalid_count"] == 1
        assert result["summary"]["valid_rate"] == pytest.approx(2 / 3)

    def test_empty_batch(self):
        result = validate_batch([])
        assert result["summary"]["total"] == 0
        assert result["summary"]["valid_rate"] == 0

    def test_all_valid(self):
        result = validate_batch(["CCO", "CC", "C"])
        assert result["summary"]["valid_count"] == 3
        assert result["summary"]["invalid_count"] == 0


class TestFilterCompounds:
    def test_filters_invalid(self):
        compounds = [
            {"name": "ethanol", "smiles": "CCO"},
            {"name": "bad", "smiles": "NOT_VALID"},
        ]
        valid, rejected = filter_valid_compounds(compounds)
        assert len(valid) == 1
        assert len(rejected) == 1
        assert valid[0]["name"] == "ethanol"
        assert rejected[0]["name"] == "bad"

    def test_adds_validation_metadata(self):
        compounds = [{"name": "ethanol", "smiles": "CCO"}]
        valid, _ = filter_valid_compounds(compounds)
        assert valid[0]["_smiles_validation"] == "passed"
