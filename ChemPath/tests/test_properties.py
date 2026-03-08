"""Tests for molecular property estimation."""

import pytest
from chempath.chemistry.properties import (
    compute_properties,
    compute_risk_scores,
    estimate_molecular_weight,
    estimate_hba,
    MolecularProperties,
)


class TestComputeProperties:
    def test_returns_properties(self):
        props = compute_properties("CCO")
        assert isinstance(props, MolecularProperties)
        assert props.estimated_mw > 0
        assert props.lipinski_violations >= 0
        assert 0.0 <= props.risk_score <= 1.0

    def test_small_molecule_low_risk(self):
        # Ethanol — very small, should have 0 Lipinski violations
        props = compute_properties("CCO")
        assert props.lipinski_violations == 0
        assert props.risk_score == 0.0

    def test_large_molecule_higher_risk(self):
        # Very long SMILES string — should flag some issues
        big_smiles = "C" * 100 + "O" * 20 + "N" * 10
        props = compute_properties(big_smiles)
        assert props.estimated_mw > 500
        assert props.lipinski_violations > 0
        assert props.risk_score > 0

    def test_risk_flags_populated(self):
        # Large molecule should have flags
        big_smiles = "C" * 100 + "O" * 20
        props = compute_properties(big_smiles)
        assert len(props.risk_flags) > 0

    def test_real_drug_smiles(self):
        # Gefitinib
        smiles = "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1"
        props = compute_properties(smiles)
        assert props.estimated_mw > 100
        assert isinstance(props.estimated_logp, float)


class TestEstimateHBA:
    def test_counts_n_and_o(self):
        assert estimate_hba("CCNCC") > 0  # has N
        assert estimate_hba("CCOCC") > 0  # has O
        assert estimate_hba("CCCC") == 0  # no N or O


class TestComputeRiskScores:
    def test_returns_dict(self):
        compounds = [
            {"chembl_id": "C1", "smiles": "CCO"},
            {"chembl_id": "C2", "smiles": "CCCCCCCCCCCCCCCCCCCC"},
        ]
        scores = compute_risk_scores(compounds)
        assert "C1" in scores
        assert "C2" in scores
        assert 0.0 <= scores["C1"] <= 1.0

    def test_missing_smiles(self):
        compounds = [{"chembl_id": "C1", "smiles": ""}]
        scores = compute_risk_scores(compounds)
        assert scores["C1"] == 0.5  # unknown
