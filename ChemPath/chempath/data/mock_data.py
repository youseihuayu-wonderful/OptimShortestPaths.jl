"""
Mock ChEMBL-style data for development and testing.
Mimics ChEMBL 36 structure with anti-cancer targets.
"""

MOCK_DATA = {
    "compounds": [
        {"chembl_id": "CHEMBL1201585", "name": "Osimertinib", "smiles": "C=CC(=O)Nc1cc(OC)c(Nc2nccc(-c3cn(C)c4ccccc34)n2)cc1N(C)CCN(C)C", "phase": 4},
        {"chembl_id": "CHEMBL1421", "name": "Gefitinib", "smiles": "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1", "phase": 4},
        {"chembl_id": "CHEMBL553", "name": "Erlotinib", "smiles": "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1", "phase": 4},
        {"chembl_id": "CHEMBL1336", "name": "Lapatinib", "smiles": "CS(=O)(=O)CCNCc1ccc(-c2ccc3ncnc(Nc4ccc(OCc5cccc(F)c5)c(Cl)c4)c3c2)o1", "phase": 4},
        {"chembl_id": "CHEMBL3545252", "name": "Alectinib", "smiles": "CCC(CC)c1cc2c(C#N)c(Nc3ccc(C(=O)c4ccncc4)cc3)nc(N)c2cn1", "phase": 4},
        {"chembl_id": "CHEMBL2007641", "name": "Ceritinib", "smiles": "CC(C)Oc1cc(C2CCNCC2)c(OC)cc1Nc1ncc(Cl)c(Nc2ccccc2S(=O)(=O)C(C)C)n1", "phase": 4},
        {"chembl_id": "CHEMBL1229517", "name": "Vemurafenib", "smiles": "CCCS(=O)(=O)Nc1ccc(-c2nc(-c3ccc(OC(F)(F)F)cc3)c(-c3ccnc(N)n3)[nH]2)cc1F", "phase": 4},
        {"chembl_id": "CHEMBL2028663", "name": "Dabrafenib", "smiles": "CC(C)(C)c1nc(-c2cccc(NS(=O)(=O)c3c(F)cccc3F)c2F)c(-c2ccnc(N)n2)s1", "phase": 4},
        {"chembl_id": "CHEMBL_INVALID", "name": "BadCompound", "smiles": "NOT_A_SMILES_XYZ", "phase": 0},
        {"chembl_id": "CHEMBL_EMPTY", "name": "EmptySmiles", "smiles": "", "phase": 0},
    ],
    "targets": [
        {"chembl_id": "CHEMBL203", "name": "EGFR", "organism": "Homo sapiens", "type": "SINGLE PROTEIN"},
        {"chembl_id": "CHEMBL1824", "name": "HER2", "organism": "Homo sapiens", "type": "SINGLE PROTEIN"},
        {"chembl_id": "CHEMBL4282", "name": "ALK", "organism": "Homo sapiens", "type": "SINGLE PROTEIN"},
        {"chembl_id": "CHEMBL5145", "name": "BRAF", "organism": "Homo sapiens", "type": "SINGLE PROTEIN"},
    ],
    "bioactivities": [
        {"compound": "CHEMBL1201585", "target": "CHEMBL203", "type": "IC50", "value": 0.5, "units": "nM", "source": "experimental"},
        {"compound": "CHEMBL1421", "target": "CHEMBL203", "type": "IC50", "value": 33.0, "units": "nM", "source": "experimental"},
        {"compound": "CHEMBL553", "target": "CHEMBL203", "type": "IC50", "value": 2.0, "units": "nM", "source": "experimental"},
        {"compound": "CHEMBL1336", "target": "CHEMBL203", "type": "IC50", "value": 10.2, "units": "nM", "source": "experimental"},
        {"compound": "CHEMBL1336", "target": "CHEMBL1824", "type": "IC50", "value": 9.8, "units": "nM", "source": "experimental"},
        {"compound": "CHEMBL1421", "target": "CHEMBL1824", "type": "IC50", "value": 500.0, "units": "nM", "source": "experimental"},
        {"compound": "CHEMBL3545252", "target": "CHEMBL4282", "type": "IC50", "value": 1.9, "units": "nM", "source": "experimental"},
        {"compound": "CHEMBL2007641", "target": "CHEMBL4282", "type": "IC50", "value": 0.2, "units": "nM", "source": "experimental"},
        {"compound": "CHEMBL1229517", "target": "CHEMBL5145", "type": "IC50", "value": 31.0, "units": "nM", "source": "experimental"},
        {"compound": "CHEMBL2028663", "target": "CHEMBL5145", "type": "IC50", "value": 0.8, "units": "nM", "source": "experimental"},
        {"compound": "CHEMBL1201585", "target": "CHEMBL1824", "type": "IC50", "value": 280.0, "units": "nM", "source": "experimental"},
        {"compound": "CHEMBL553", "target": "CHEMBL1824", "type": "IC50", "value": 350.0, "units": "nM", "source": "experimental"},
    ],
    "toxicity": {
        "CHEMBL1201585": {"hepatotoxicity": 0.15, "cardiotoxicity": 0.20, "overall": 0.18},
        "CHEMBL1421": {"hepatotoxicity": 0.25, "cardiotoxicity": 0.10, "overall": 0.20},
        "CHEMBL553": {"hepatotoxicity": 0.20, "cardiotoxicity": 0.15, "overall": 0.18},
        "CHEMBL1336": {"hepatotoxicity": 0.30, "cardiotoxicity": 0.35, "overall": 0.33},
        "CHEMBL3545252": {"hepatotoxicity": 0.10, "cardiotoxicity": 0.12, "overall": 0.11},
        "CHEMBL2007641": {"hepatotoxicity": 0.22, "cardiotoxicity": 0.18, "overall": 0.20},
        "CHEMBL1229517": {"hepatotoxicity": 0.40, "cardiotoxicity": 0.25, "overall": 0.35},
        "CHEMBL2028663": {"hepatotoxicity": 0.15, "cardiotoxicity": 0.20, "overall": 0.18},
    },
}


def load_mock_data() -> dict:
    """Load mock ChEMBL-style data."""
    return MOCK_DATA
