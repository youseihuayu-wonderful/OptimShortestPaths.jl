"""
Curated biological network data for multi-hop shortest-path optimization.

This encodes well-established biological relationships:
  - Signaling pathways (KEGG/Reactome-derived)
  - Pathway-disease associations (DisGeNET/KEGG Disease-derived)
  - Protein-protein interactions (STRING/BioGRID-derived)

All relationships are curated from literature; weights reflect
strength/confidence of the biological association.

Weight semantics (lower = stronger connection):
  - PPI confidence 0.9 → weight 0.1  (strong interaction = low cost to traverse)
  - Pathway centrality "high" → weight 0.5  (key driver = easy to reach)
  - Disease association "strong" → weight 0.3  (well-validated link)
"""

# ---------------------------------------------------------------------------
# Signaling pathways relevant to anti-cancer targets
# ---------------------------------------------------------------------------

PATHWAYS = [
    {
        "id": "PATH_ERBB",
        "name": "ErbB Signaling",
        "description": "EGFR/HER family receptor tyrosine kinase signaling",
        "kegg_id": "hsa04012",
    },
    {
        "id": "PATH_MAPK",
        "name": "MAPK/ERK Cascade",
        "description": "RAS-RAF-MEK-ERK mitogen-activated protein kinase cascade",
        "kegg_id": "hsa04010",
    },
    {
        "id": "PATH_PI3K_AKT",
        "name": "PI3K/AKT/mTOR",
        "description": "Phosphoinositide 3-kinase / AKT / mTOR survival signaling",
        "kegg_id": "hsa04151",
    },
    {
        "id": "PATH_ANGIOGENESIS",
        "name": "VEGF/Angiogenesis",
        "description": "Vascular endothelial growth factor angiogenesis pathway",
        "kegg_id": "hsa04370",
    },
    {
        "id": "PATH_HGF_MET",
        "name": "HGF/MET Signaling",
        "description": "Hepatocyte growth factor / MET receptor invasion signaling",
        "kegg_id": "hsa04360",
    },
    {
        "id": "PATH_FGF",
        "name": "FGF Signaling",
        "description": "Fibroblast growth factor receptor signaling",
        "kegg_id": "hsa04010",
    },
    {
        "id": "PATH_BCR_ABL",
        "name": "BCR-ABL Signaling",
        "description": "BCR-ABL fusion oncoprotein signaling in CML",
        "kegg_id": "hsa05220",
    },
    {
        "id": "PATH_APOPTOSIS",
        "name": "Apoptosis Regulation",
        "description": "Programmed cell death / survival balance",
        "kegg_id": "hsa04210",
    },
    {
        "id": "PATH_CELL_CYCLE",
        "name": "Cell Cycle Control",
        "description": "Cell cycle progression and checkpoint regulation",
        "kegg_id": "hsa04110",
    },
    {
        "id": "PATH_EMT",
        "name": "EMT / Cell Migration",
        "description": "Epithelial-mesenchymal transition and cell motility",
        "kegg_id": "hsa04310",
    },
]

# ---------------------------------------------------------------------------
# Diseases (cancers) linked to these pathways
# ---------------------------------------------------------------------------

DISEASES = [
    {
        "id": "DIS_NSCLC",
        "name": "NSCLC",
        "full_name": "Non-Small Cell Lung Cancer",
        "icd10": "C34",
    },
    {
        "id": "DIS_BREAST",
        "name": "Breast Cancer",
        "full_name": "Breast Carcinoma (HER2+)",
        "icd10": "C50",
    },
    {
        "id": "DIS_CML",
        "name": "CML",
        "full_name": "Chronic Myeloid Leukemia",
        "icd10": "C92.1",
    },
    {
        "id": "DIS_MELANOMA",
        "name": "Melanoma",
        "full_name": "Malignant Melanoma",
        "icd10": "C43",
    },
    {
        "id": "DIS_CRC",
        "name": "Colorectal Cancer",
        "full_name": "Colorectal Carcinoma",
        "icd10": "C18",
    },
    {
        "id": "DIS_RCC",
        "name": "Renal Cell Carcinoma",
        "full_name": "Renal Cell Carcinoma",
        "icd10": "C64",
    },
    {
        "id": "DIS_GBM",
        "name": "Glioblastoma",
        "full_name": "Glioblastoma Multiforme",
        "icd10": "C71",
    },
    {
        "id": "DIS_HCC",
        "name": "Hepatocellular Carcinoma",
        "full_name": "Hepatocellular Carcinoma",
        "icd10": "C22.0",
    },
    {
        "id": "DIS_GIST",
        "name": "GIST",
        "full_name": "Gastrointestinal Stromal Tumor",
        "icd10": "C49",
    },
    {
        "id": "DIS_THYROID",
        "name": "Thyroid Cancer",
        "full_name": "Differentiated Thyroid Carcinoma",
        "icd10": "C73",
    },
]

# ---------------------------------------------------------------------------
# Target → Pathway edges
# Weight = inverse of centrality (lower = target is more central to pathway)
# A target can participate in multiple pathways
# ---------------------------------------------------------------------------

TARGET_PATHWAY_EDGES = [
    # EGFR (CHEMBL203) — central to ErbB, feeds into MAPK and PI3K
    {"source": "CHEMBL203", "target": "PATH_ERBB", "weight": 0.3, "role": "receptor", "evidence": "KEGG hsa04012"},
    {"source": "CHEMBL203", "target": "PATH_MAPK", "weight": 0.8, "role": "upstream_activator", "evidence": "EGFR→RAS→RAF→MEK→ERK"},
    {"source": "CHEMBL203", "target": "PATH_PI3K_AKT", "weight": 0.9, "role": "upstream_activator", "evidence": "EGFR→PI3K→AKT"},

    # HER2 (CHEMBL1824) — ErbB family, feeds into PI3K primarily
    {"source": "CHEMBL1824", "target": "PATH_ERBB", "weight": 0.3, "role": "receptor", "evidence": "KEGG hsa04012"},
    {"source": "CHEMBL1824", "target": "PATH_PI3K_AKT", "weight": 0.7, "role": "upstream_activator", "evidence": "HER2→PI3K preferred"},
    {"source": "CHEMBL1824", "target": "PATH_MAPK", "weight": 0.9, "role": "upstream_activator", "evidence": "HER2→SHC→RAS"},

    # ALK (CHEMBL4282) — feeds into MAPK, PI3K, and cell cycle
    {"source": "CHEMBL4282", "target": "PATH_MAPK", "weight": 0.5, "role": "upstream_activator", "evidence": "ALK→RAS→MAPK"},
    {"source": "CHEMBL4282", "target": "PATH_PI3K_AKT", "weight": 0.6, "role": "upstream_activator", "evidence": "ALK→PI3K→AKT"},
    {"source": "CHEMBL4282", "target": "PATH_CELL_CYCLE", "weight": 0.8, "role": "indirect_regulator", "evidence": "ALK→STAT3→Cyclin D"},

    # BRAF (CHEMBL5145) — central node in MAPK cascade
    {"source": "CHEMBL5145", "target": "PATH_MAPK", "weight": 0.2, "role": "central_kinase", "evidence": "RAF is the core of MAPK cascade"},
    {"source": "CHEMBL5145", "target": "PATH_CELL_CYCLE", "weight": 0.7, "role": "indirect_regulator", "evidence": "MAPK→Cyclin D1"},

    # ABL1 (CHEMBL1862) — BCR-ABL signaling, also MAPK
    {"source": "CHEMBL1862", "target": "PATH_BCR_ABL", "weight": 0.1, "role": "oncogenic_driver", "evidence": "BCR-ABL fusion = constitutive activation"},
    {"source": "CHEMBL1862", "target": "PATH_MAPK", "weight": 0.7, "role": "upstream_activator", "evidence": "ABL→RAS→MAPK"},
    {"source": "CHEMBL1862", "target": "PATH_PI3K_AKT", "weight": 0.8, "role": "upstream_activator", "evidence": "ABL→PI3K"},

    # VEGFR2 (CHEMBL2842) — angiogenesis, PI3K cross-talk
    {"source": "CHEMBL2842", "target": "PATH_ANGIOGENESIS", "weight": 0.2, "role": "central_receptor", "evidence": "VEGFR2 = primary VEGF receptor"},
    {"source": "CHEMBL2842", "target": "PATH_PI3K_AKT", "weight": 0.8, "role": "upstream_activator", "evidence": "VEGFR2→PI3K→AKT→eNOS"},
    {"source": "CHEMBL2842", "target": "PATH_MAPK", "weight": 0.9, "role": "upstream_activator", "evidence": "VEGFR2→PLCγ→PKC→MAPK"},

    # MET (CHEMBL4005) — HGF/MET, invasion, EMT
    {"source": "CHEMBL4005", "target": "PATH_HGF_MET", "weight": 0.2, "role": "receptor", "evidence": "MET = HGF receptor"},
    {"source": "CHEMBL4005", "target": "PATH_EMT", "weight": 0.5, "role": "driver", "evidence": "MET promotes EMT and invasion"},
    {"source": "CHEMBL4005", "target": "PATH_MAPK", "weight": 0.7, "role": "upstream_activator", "evidence": "MET→GRB2→SOS→RAS"},
    {"source": "CHEMBL4005", "target": "PATH_PI3K_AKT", "weight": 0.7, "role": "upstream_activator", "evidence": "MET→PI3K→AKT"},

    # FGFR1 (CHEMBL3594) — FGF signaling, MAPK
    {"source": "CHEMBL3594", "target": "PATH_FGF", "weight": 0.2, "role": "receptor", "evidence": "FGFR1 = FGF receptor"},
    {"source": "CHEMBL3594", "target": "PATH_MAPK", "weight": 0.6, "role": "upstream_activator", "evidence": "FGFR→FRS2→RAS→MAPK"},
    {"source": "CHEMBL3594", "target": "PATH_PI3K_AKT", "weight": 0.8, "role": "upstream_activator", "evidence": "FGFR→PI3K"},

    # PI3Kα (CHEMBL4630) — central to PI3K/AKT/mTOR
    {"source": "CHEMBL4630", "target": "PATH_PI3K_AKT", "weight": 0.1, "role": "central_kinase", "evidence": "PI3Kα is the catalytic subunit"},
    {"source": "CHEMBL4630", "target": "PATH_APOPTOSIS", "weight": 0.6, "role": "survival_signal", "evidence": "PI3K→AKT→BAD (anti-apoptotic)"},

    # mTOR (CHEMBL2185) — PI3K/AKT/mTOR, cell growth
    {"source": "CHEMBL2185", "target": "PATH_PI3K_AKT", "weight": 0.3, "role": "downstream_effector", "evidence": "AKT→mTOR→S6K/4EBP1"},
    {"source": "CHEMBL2185", "target": "PATH_CELL_CYCLE", "weight": 0.5, "role": "growth_regulator", "evidence": "mTOR→S6K→cell growth"},
    {"source": "CHEMBL2185", "target": "PATH_APOPTOSIS", "weight": 0.7, "role": "survival_signal", "evidence": "mTOR→autophagy regulation"},
]

# ---------------------------------------------------------------------------
# Pathway → Disease edges
# Weight = inverse of association strength (lower = stronger link)
# ---------------------------------------------------------------------------

PATHWAY_DISEASE_EDGES = [
    # ErbB signaling → cancers
    {"source": "PATH_ERBB", "target": "DIS_NSCLC", "weight": 0.3, "evidence": "EGFR mutations in 15-30% of NSCLC"},
    {"source": "PATH_ERBB", "target": "DIS_BREAST", "weight": 0.3, "evidence": "HER2 amplified in 20% of breast cancer"},
    {"source": "PATH_ERBB", "target": "DIS_CRC", "weight": 0.6, "evidence": "EGFR overexpression in CRC"},
    {"source": "PATH_ERBB", "target": "DIS_GBM", "weight": 0.5, "evidence": "EGFRvIII in ~30% of GBM"},

    # MAPK cascade → cancers
    {"source": "PATH_MAPK", "target": "DIS_MELANOMA", "weight": 0.2, "evidence": "BRAF V600E in 50% of melanoma"},
    {"source": "PATH_MAPK", "target": "DIS_NSCLC", "weight": 0.5, "evidence": "KRAS mutations in 25% of NSCLC"},
    {"source": "PATH_MAPK", "target": "DIS_CRC", "weight": 0.4, "evidence": "KRAS/BRAF mutations in 50% of CRC"},
    {"source": "PATH_MAPK", "target": "DIS_THYROID", "weight": 0.4, "evidence": "BRAF V600E in 45% of PTC"},

    # PI3K/AKT/mTOR → cancers
    {"source": "PATH_PI3K_AKT", "target": "DIS_BREAST", "weight": 0.3, "evidence": "PIK3CA mutations in 35% of HR+ breast cancer"},
    {"source": "PATH_PI3K_AKT", "target": "DIS_RCC", "weight": 0.5, "evidence": "mTOR pathway activation in RCC"},
    {"source": "PATH_PI3K_AKT", "target": "DIS_GBM", "weight": 0.5, "evidence": "PTEN loss / PI3K activation in GBM"},
    {"source": "PATH_PI3K_AKT", "target": "DIS_HCC", "weight": 0.6, "evidence": "PI3K/AKT/mTOR in HCC progression"},
    {"source": "PATH_PI3K_AKT", "target": "DIS_CRC", "weight": 0.6, "evidence": "PIK3CA mutations in CRC"},

    # Angiogenesis → cancers (angiogenesis is pan-cancer but especially relevant for these)
    {"source": "PATH_ANGIOGENESIS", "target": "DIS_RCC", "weight": 0.2, "evidence": "VHL→VEGF is the primary driver of clear cell RCC"},
    {"source": "PATH_ANGIOGENESIS", "target": "DIS_HCC", "weight": 0.4, "evidence": "Hypervascular tumor, VEGF-driven"},
    {"source": "PATH_ANGIOGENESIS", "target": "DIS_CRC", "weight": 0.5, "evidence": "Anti-VEGF (bevacizumab) approved for mCRC"},
    {"source": "PATH_ANGIOGENESIS", "target": "DIS_GBM", "weight": 0.5, "evidence": "Highly vascular tumor"},

    # HGF/MET → cancers
    {"source": "PATH_HGF_MET", "target": "DIS_NSCLC", "weight": 0.4, "evidence": "MET amplification in EGFR-resistant NSCLC"},
    {"source": "PATH_HGF_MET", "target": "DIS_HCC", "weight": 0.3, "evidence": "MET overexpression in HCC, HGF from liver stroma"},
    {"source": "PATH_HGF_MET", "target": "DIS_RCC", "weight": 0.5, "evidence": "MET/HGF in papillary RCC"},
    {"source": "PATH_HGF_MET", "target": "DIS_GBM", "weight": 0.6, "evidence": "MET amplification in GBM"},

    # FGF signaling → cancers
    {"source": "PATH_FGF", "target": "DIS_BREAST", "weight": 0.5, "evidence": "FGFR1 amplification in ER+ breast cancer"},
    {"source": "PATH_FGF", "target": "DIS_NSCLC", "weight": 0.6, "evidence": "FGFR alterations in squamous NSCLC"},
    {"source": "PATH_FGF", "target": "DIS_HCC", "weight": 0.7, "evidence": "FGF19 amplification in HCC subset"},

    # BCR-ABL → CML (very strong, nearly 1:1)
    {"source": "PATH_BCR_ABL", "target": "DIS_CML", "weight": 0.1, "evidence": "BCR-ABL fusion defines CML (Philadelphia chromosome)"},
    {"source": "PATH_BCR_ABL", "target": "DIS_GIST", "weight": 0.7, "evidence": "ABL-related kinase signaling in some GIST"},

    # Apoptosis → pan-cancer
    {"source": "PATH_APOPTOSIS", "target": "DIS_CML", "weight": 0.5, "evidence": "BCR-ABL blocks apoptosis"},
    {"source": "PATH_APOPTOSIS", "target": "DIS_BREAST", "weight": 0.6, "evidence": "Apoptosis evasion in breast cancer"},

    # Cell cycle → cancers
    {"source": "PATH_CELL_CYCLE", "target": "DIS_MELANOMA", "weight": 0.5, "evidence": "CDKN2A loss in melanoma"},
    {"source": "PATH_CELL_CYCLE", "target": "DIS_NSCLC", "weight": 0.6, "evidence": "Cell cycle dysregulation in NSCLC"},
    {"source": "PATH_CELL_CYCLE", "target": "DIS_BREAST", "weight": 0.5, "evidence": "CDK4/6 pathway in HR+ breast cancer"},

    # EMT → metastatic cancers
    {"source": "PATH_EMT", "target": "DIS_HCC", "weight": 0.4, "evidence": "EMT drives HCC invasion and metastasis"},
    {"source": "PATH_EMT", "target": "DIS_NSCLC", "weight": 0.5, "evidence": "EMT in NSCLC metastasis"},
    {"source": "PATH_EMT", "target": "DIS_BREAST", "weight": 0.5, "evidence": "EMT in triple-negative breast cancer"},
]

# ---------------------------------------------------------------------------
# Protein-protein interactions (PPIs) between our targets
# Weight = 1 - confidence (STRING-derived; higher confidence = lower weight)
# These are directed for signaling cascades (A activates/inhibits B)
# ---------------------------------------------------------------------------

PPI_EDGES = [
    # EGFR ↔ HER2 heterodimerization
    {"source": "CHEMBL203", "target": "CHEMBL1824", "weight": 0.2, "type": "heterodimerization",
     "evidence": "EGFR/HER2 heterodimer is the most potent ErbB dimer", "bidirectional": True},

    # EGFR → MET cross-talk (MET amplification = EGFR resistance)
    {"source": "CHEMBL203", "target": "CHEMBL4005", "weight": 0.5, "type": "cross_resistance",
     "evidence": "MET amplification bypasses EGFR inhibition", "bidirectional": False},

    # BRAF → MAPK is inherently captured in target→pathway, but BRAF→MEK is a direct PPI
    # BRAF and ABL1 are both kinases that converge on MAPK
    {"source": "CHEMBL5145", "target": "CHEMBL1862", "weight": 0.8, "type": "pathway_convergence",
     "evidence": "Both converge on MAPK but through different branches", "bidirectional": False},

    # VEGFR2 → PI3Kα direct activation
    {"source": "CHEMBL2842", "target": "CHEMBL4630", "weight": 0.4, "type": "direct_activation",
     "evidence": "VEGFR2 directly phosphorylates PI3K", "bidirectional": False},

    # MET → PI3Kα activation
    {"source": "CHEMBL4005", "target": "CHEMBL4630", "weight": 0.5, "type": "direct_activation",
     "evidence": "MET activates PI3K through GAB1 adaptor", "bidirectional": False},

    # PI3Kα → mTOR (PI3K/AKT/mTOR cascade)
    {"source": "CHEMBL4630", "target": "CHEMBL2185", "weight": 0.3, "type": "signaling_cascade",
     "evidence": "PI3K→AKT→TSC1/2→mTOR activation", "bidirectional": False},

    # FGFR1 → PI3Kα cross-talk
    {"source": "CHEMBL3594", "target": "CHEMBL4630", "weight": 0.6, "type": "direct_activation",
     "evidence": "FGFR activates PI3K through FRS2/GAB1", "bidirectional": False},

    # HER2 → PI3Kα (strong direct interaction)
    {"source": "CHEMBL1824", "target": "CHEMBL4630", "weight": 0.3, "type": "direct_activation",
     "evidence": "HER2/HER3 dimer directly activates PI3K", "bidirectional": False},

    # ALK → ABL1 (both in overlapping oncogenic programs)
    {"source": "CHEMBL4282", "target": "CHEMBL1862", "weight": 0.9, "type": "functional_overlap",
     "evidence": "ALK and ABL share downstream substrates (STAT3, CRK)", "bidirectional": False},

    # EGFR → FGFR1 (resistance mechanism)
    {"source": "CHEMBL203", "target": "CHEMBL3594", "weight": 0.6, "type": "cross_resistance",
     "evidence": "FGFR activation as EGFR TKI resistance mechanism", "bidirectional": False},

    # mTOR feedback → PI3Kα (negative feedback loop)
    {"source": "CHEMBL2185", "target": "CHEMBL4630", "weight": 0.5, "type": "negative_feedback",
     "evidence": "mTOR/S6K→IRS1 degradation reduces PI3K input (paradoxical activation on mTOR inhibition)",
     "bidirectional": False},
]


def get_biological_network() -> dict:
    """Return the complete curated biological network data."""
    return {
        "pathways": PATHWAYS,
        "diseases": DISEASES,
        "target_pathway_edges": TARGET_PATHWAY_EDGES,
        "pathway_disease_edges": PATHWAY_DISEASE_EDGES,
        "ppi_edges": PPI_EDGES,
    }
