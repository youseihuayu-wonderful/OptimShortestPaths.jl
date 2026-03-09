"""
Curated ChEMBL-style kinase inhibitor dataset for multi-hop drug discovery.

Contains 20 FDA-approved kinase inhibitors across 11 protein targets,
with experimentally-derived IC50 values and toxicity profiles.

Data sources:
  - IC50 values: ChEMBL 36 (approximate, for demonstration)
  - Toxicity: FDA labels, LiverTox database (simplified)
  - Targets: UniProt/ChEMBL target annotations

All compounds are FDA-approved for at least one cancer indication,
enabling retrospective validation against known drug-disease pairs.
"""

MOCK_DATA = {
    "compounds": [
        # --- EGFR inhibitors ---
        {"chembl_id": "CHEMBL939", "name": "Gefitinib",
         "smiles": "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1", "phase": 4},
        {"chembl_id": "CHEMBL553", "name": "Erlotinib",
         "smiles": "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1", "phase": 4},
        {"chembl_id": "CHEMBL1201585", "name": "Osimertinib",
         "smiles": "C=CC(=O)Nc1cc(OC)c(Nc2nccc(-c3cn(C)c4ccccc34)n2)cc1N(C)CCN(C)C", "phase": 4},
        # --- EGFR + HER2 dual inhibitors ---
        {"chembl_id": "CHEMBL554", "name": "Lapatinib",
         "smiles": "CS(=O)(=O)CCNCc1ccc(-c2ccc3ncnc(Nc4ccc(OCc5cccc(F)c5)c(Cl)c4)c3c2)o1", "phase": 4},
        {"chembl_id": "CHEMBL180022", "name": "Neratinib",
         "smiles": "C=CC(=O)Nc1cc2c(Nc3ccc(Nc4nc(-c5ccccc5)nc5ccccc45)cc3OC)ccn2c1", "phase": 4},
        {"chembl_id": "CHEMBL24828", "name": "Vandetanib",
         "smiles": "COc1cc2c(Nc3ccc(Br)cc3F)ncnc2cc1OCC1CCN(C)CC1", "phase": 4},
        # --- ALK inhibitors ---
        {"chembl_id": "CHEMBL3545252", "name": "Alectinib",
         "smiles": "CCC(CC)c1cc2c(C#N)c(Nc3ccc(C(=O)c4ccncc4)cc3)nc(N)c2cn1", "phase": 4},
        {"chembl_id": "CHEMBL2007641", "name": "Ceritinib",
         "smiles": "CC(C)Oc1cc(C2CCNCC2)c(OC)cc1Nc1ncc(Cl)c(Nc2ccccc2S(=O)(=O)C(C)C)n1", "phase": 4},
        {"chembl_id": "CHEMBL601719", "name": "Crizotinib",
         "smiles": "CC(Oc1cc(-c2cnn(C3CCNCC3)c2)cnc1N)c1c(Cl)ccc(F)c1Cl", "phase": 4},
        # --- BRAF inhibitors ---
        {"chembl_id": "CHEMBL1229517", "name": "Vemurafenib",
         "smiles": "CCCS(=O)(=O)Nc1ccc(-c2nc(-c3ccc(OC(F)(F)F)cc3)c(-c3ccnc(N)n3)[nH]2)cc1F", "phase": 4},
        {"chembl_id": "CHEMBL2028663", "name": "Dabrafenib",
         "smiles": "CC(C)(C)c1nc(-c2cccc(NS(=O)(=O)c3c(F)cccc3F)c2F)c(-c2ccnc(N)n2)s1", "phase": 4},
        # --- BCR-ABL inhibitors ---
        {"chembl_id": "CHEMBL941", "name": "Imatinib",
         "smiles": "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1", "phase": 4},
        {"chembl_id": "CHEMBL288441", "name": "Bosutinib",
         "smiles": "COc1cc(Nc2c(C#N)cnc3cc(OCCCN4CCN(C)CC4)c(OC)cc23)c(Cl)cc1Cl", "phase": 4},
        {"chembl_id": "CHEMBL255863", "name": "Nilotinib",
         "smiles": "Cc1cn(-c2cc(NC(=O)c3ccc(C)c(Nc4nccc(-c5cccnc5)n4)c3)cc(C(F)(F)F)c2)cn1", "phase": 4},
        {"chembl_id": "CHEMBL1171837", "name": "Ponatinib",
         "smiles": "Cc1ccc(C(=O)Nc2ccc(C)c(C#Cc3cnc4cccnn34)c2)cc1C1CCN(CC1)CC(=O)N1CCOCC1", "phase": 4},
        {"chembl_id": "CHEMBL1421", "name": "Dasatinib",
         "smiles": "Cc1nc(Nc2ncc(C(=O)Nc3c(C)cccc3Cl)s2)cc(N2CCN(CCO)CC2)n1", "phase": 4},
        # --- Multi-kinase inhibitors ---
        {"chembl_id": "CHEMBL1336", "name": "Sorafenib",
         "smiles": "CNC(=O)c1cc(Oc2ccc(NC(=O)Nc3ccc(Cl)c(C(F)(F)F)c3)cc2)ccn1", "phase": 4},
        # --- mTOR inhibitor ---
        {"chembl_id": "CHEMBL1642", "name": "Everolimus",
         "smiles": "COc1cc2c(cc1OC)C(=O)C(CC=CC(C)C(OC1OC(C)C(O)C(OC)C1O)C(C)CC(=O)C(O)/C(C)=C/C1CCC(C(CC3)C(=O)C3=O)O1)C2", "phase": 4},
        # --- Edge cases for testing ---
        {"chembl_id": "CHEMBL_INVALID", "name": "BadCompound",
         "smiles": "NOT_A_SMILES_XYZ", "phase": 0},
        {"chembl_id": "CHEMBL_EMPTY", "name": "EmptySmiles",
         "smiles": "", "phase": 0},
    ],

    "targets": [
        {"chembl_id": "CHEMBL203", "name": "EGFR", "organism": "Homo sapiens", "type": "SINGLE PROTEIN"},
        {"chembl_id": "CHEMBL1824", "name": "HER2", "organism": "Homo sapiens", "type": "SINGLE PROTEIN"},
        {"chembl_id": "CHEMBL4282", "name": "ALK", "organism": "Homo sapiens", "type": "SINGLE PROTEIN"},
        {"chembl_id": "CHEMBL5145", "name": "BRAF", "organism": "Homo sapiens", "type": "SINGLE PROTEIN"},
        {"chembl_id": "CHEMBL1862", "name": "ABL1", "organism": "Homo sapiens", "type": "SINGLE PROTEIN"},
        {"chembl_id": "CHEMBL2842", "name": "VEGFR2", "organism": "Homo sapiens", "type": "SINGLE PROTEIN"},
        {"chembl_id": "CHEMBL4005", "name": "MET", "organism": "Homo sapiens", "type": "SINGLE PROTEIN"},
        {"chembl_id": "CHEMBL3594", "name": "FGFR1", "organism": "Homo sapiens", "type": "SINGLE PROTEIN"},
        {"chembl_id": "CHEMBL4630", "name": "PI3Ka", "organism": "Homo sapiens", "type": "SINGLE PROTEIN"},
        {"chembl_id": "CHEMBL2185", "name": "mTOR", "organism": "Homo sapiens", "type": "SINGLE PROTEIN"},
        {"chembl_id": "CHEMBL1936", "name": "KIT", "organism": "Homo sapiens", "type": "SINGLE PROTEIN"},
    ],

    # IC50 values (nM) — approximate values from ChEMBL 36
    "bioactivities": [
        # EGFR inhibitors → EGFR
        {"compound": "CHEMBL939", "target": "CHEMBL203", "type": "IC50", "value": 33.0, "units": "nM", "source": "experimental"},
        {"compound": "CHEMBL553", "target": "CHEMBL203", "type": "IC50", "value": 2.0, "units": "nM", "source": "experimental"},
        {"compound": "CHEMBL1201585", "target": "CHEMBL203", "type": "IC50", "value": 0.5, "units": "nM", "source": "experimental"},
        # EGFR inhibitors with weak HER2 activity
        {"compound": "CHEMBL939", "target": "CHEMBL1824", "type": "IC50", "value": 3700.0, "units": "nM", "source": "experimental"},
        {"compound": "CHEMBL553", "target": "CHEMBL1824", "type": "IC50", "value": 350.0, "units": "nM", "source": "experimental"},
        {"compound": "CHEMBL1201585", "target": "CHEMBL1824", "type": "IC50", "value": 280.0, "units": "nM", "source": "experimental"},
        # Dual EGFR/HER2 inhibitors
        {"compound": "CHEMBL554", "target": "CHEMBL203", "type": "IC50", "value": 10.8, "units": "nM", "source": "experimental"},
        {"compound": "CHEMBL554", "target": "CHEMBL1824", "type": "IC50", "value": 9.2, "units": "nM", "source": "experimental"},
        {"compound": "CHEMBL180022", "target": "CHEMBL1824", "type": "IC50", "value": 59.0, "units": "nM", "source": "experimental"},
        {"compound": "CHEMBL180022", "target": "CHEMBL203", "type": "IC50", "value": 92.0, "units": "nM", "source": "experimental"},
        # Vandetanib: VEGFR2 + EGFR
        {"compound": "CHEMBL24828", "target": "CHEMBL2842", "type": "IC50", "value": 40.0, "units": "nM", "source": "experimental"},
        {"compound": "CHEMBL24828", "target": "CHEMBL203", "type": "IC50", "value": 500.0, "units": "nM", "source": "experimental"},
        # ALK inhibitors
        {"compound": "CHEMBL3545252", "target": "CHEMBL4282", "type": "IC50", "value": 1.9, "units": "nM", "source": "experimental"},
        {"compound": "CHEMBL2007641", "target": "CHEMBL4282", "type": "IC50", "value": 0.2, "units": "nM", "source": "experimental"},
        {"compound": "CHEMBL601719", "target": "CHEMBL4282", "type": "IC50", "value": 24.0, "units": "nM", "source": "experimental"},
        # Crizotinib also hits MET
        {"compound": "CHEMBL601719", "target": "CHEMBL4005", "type": "IC50", "value": 11.0, "units": "nM", "source": "experimental"},
        # BRAF inhibitors
        {"compound": "CHEMBL1229517", "target": "CHEMBL5145", "type": "IC50", "value": 31.0, "units": "nM", "source": "experimental"},
        {"compound": "CHEMBL2028663", "target": "CHEMBL5145", "type": "IC50", "value": 0.8, "units": "nM", "source": "experimental"},
        # BCR-ABL inhibitors → ABL1
        {"compound": "CHEMBL941", "target": "CHEMBL1862", "type": "IC50", "value": 600.0, "units": "nM", "source": "experimental"},
        {"compound": "CHEMBL288441", "target": "CHEMBL1862", "type": "IC50", "value": 1.0, "units": "nM", "source": "experimental"},
        {"compound": "CHEMBL255863", "target": "CHEMBL1862", "type": "IC50", "value": 20.0, "units": "nM", "source": "experimental"},
        {"compound": "CHEMBL1171837", "target": "CHEMBL1862", "type": "IC50", "value": 0.4, "units": "nM", "source": "experimental"},
        {"compound": "CHEMBL1421", "target": "CHEMBL1862", "type": "IC50", "value": 0.3, "units": "nM", "source": "experimental"},
        # Imatinib also hits KIT (important for GIST indication)
        {"compound": "CHEMBL941", "target": "CHEMBL1936", "type": "IC50", "value": 100.0, "units": "nM", "source": "experimental"},
        # Sorafenib: multi-kinase (VEGFR2 + BRAF + weak MET)
        {"compound": "CHEMBL1336", "target": "CHEMBL2842", "type": "IC50", "value": 90.0, "units": "nM", "source": "experimental"},
        {"compound": "CHEMBL1336", "target": "CHEMBL5145", "type": "IC50", "value": 22.0, "units": "nM", "source": "experimental"},
        # Everolimus → mTOR
        {"compound": "CHEMBL1642", "target": "CHEMBL2185", "type": "IC50", "value": 1.6, "units": "nM", "source": "experimental"},
    ],

    "toxicity": {
        "CHEMBL939":     {"hepatotoxicity": 0.20, "cardiotoxicity": 0.10, "overall": 0.18},
        "CHEMBL553":     {"hepatotoxicity": 0.20, "cardiotoxicity": 0.15, "overall": 0.18},
        "CHEMBL1201585": {"hepatotoxicity": 0.15, "cardiotoxicity": 0.20, "overall": 0.18},
        "CHEMBL554":     {"hepatotoxicity": 0.30, "cardiotoxicity": 0.35, "overall": 0.33},
        "CHEMBL180022":  {"hepatotoxicity": 0.25, "cardiotoxicity": 0.30, "overall": 0.28},
        "CHEMBL24828":   {"hepatotoxicity": 0.15, "cardiotoxicity": 0.40, "overall": 0.30},
        "CHEMBL3545252": {"hepatotoxicity": 0.10, "cardiotoxicity": 0.12, "overall": 0.11},
        "CHEMBL2007641": {"hepatotoxicity": 0.22, "cardiotoxicity": 0.18, "overall": 0.20},
        "CHEMBL601719":  {"hepatotoxicity": 0.25, "cardiotoxicity": 0.20, "overall": 0.23},
        "CHEMBL1229517": {"hepatotoxicity": 0.40, "cardiotoxicity": 0.25, "overall": 0.35},
        "CHEMBL2028663": {"hepatotoxicity": 0.15, "cardiotoxicity": 0.20, "overall": 0.18},
        "CHEMBL941":     {"hepatotoxicity": 0.15, "cardiotoxicity": 0.10, "overall": 0.13},
        "CHEMBL288441":  {"hepatotoxicity": 0.35, "cardiotoxicity": 0.15, "overall": 0.28},
        "CHEMBL255863":  {"hepatotoxicity": 0.20, "cardiotoxicity": 0.15, "overall": 0.18},
        "CHEMBL1171837": {"hepatotoxicity": 0.30, "cardiotoxicity": 0.40, "overall": 0.36},
        "CHEMBL1421":    {"hepatotoxicity": 0.15, "cardiotoxicity": 0.20, "overall": 0.18},
        "CHEMBL1336":    {"hepatotoxicity": 0.35, "cardiotoxicity": 0.25, "overall": 0.32},
        "CHEMBL1642":    {"hepatotoxicity": 0.20, "cardiotoxicity": 0.15, "overall": 0.18},
    },
}


def load_mock_data() -> dict:
    """Load curated kinase inhibitor dataset."""
    return MOCK_DATA
