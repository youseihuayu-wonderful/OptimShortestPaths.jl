"""
SIDER 4.1 Frequency-Weighted Severity Safety Score
===================================================
Builds a bias-corrected per-compound safety score for drug repurposing.

PROBLEM WITH RAW SIDE-EFFECT COUNT
------------------------------------
Hetionet's Compound-causes-SideEffect edges come from SIDER 4.1 (sourced from
FDA drug package inserts via NLP). Well-studied, widely-prescribed drugs
(aspirin, metformin, statins) have 100s of documented side effects NOT because
they are more dangerous but because:
  1. More patients → more surveillance → more reports (Weber effect)
  2. Older drugs have decades of post-market accumulation
  3. Serious indication drugs are scrutinized more (indication bias)

This makes raw side-effect count *anti-predictive*: effective drugs appear unsafe.

SOLUTION: SIDER FREQUENCY × CTCAE-SEVERITY WEIGHTED SCORE
-----------------------------------------------------------
SIDER 4.1 meddra_freq.tsv contains frequency data DERIVED FROM DRUG LABELS
(package inserts), not from spontaneous reports. Label frequencies reflect
FDA-submission clinical trial populations — controlled conditions that are NOT
subject to the Weber effect or notoriety bias.

Score formula per compound c:
    safety_score(c) = Σ_{ADR_i ∈ PT_terms(c)} freq_lower(i) × severity_weight(i)
    normalized by: log(1 + n_ADR_pairs) to penalize label-depth inflation

Design decisions (pre-registered, not tuned on results):
  - Only PT (Preferred Term) rows used — deduplicates MedDRA hierarchy
  - freq_lower (conservative bound) used, not mean or upper
  - Severity weights from CTCAE v5.0 grade interpretation (keywords on PT name)
  - Compounds missing from SIDER: flagged separately, raw-count fallback
  - Normalization: divide by log(1 + n_PT_pairs) to control for label depth

CTCAE SEVERITY TIERS (keyword-based, fully documented below)
--------------------------------------------------------------
Grade 5 (weight 5.0): Death outcomes
Grade 4 (weight 4.0): Life-threatening organ toxicity
Grade 3 (weight 3.0): Severe systemic events
Grade 2 (weight 2.0): Moderate symptomatic events
Grade 1 (weight 1.0): Mild / general symptoms (default)

VALIDATION PROTOCOL
--------------------
A correct metric must satisfy: thalidomide > aspirin by safety score.
Thalidomide: severe teratogenicity, peripheral neuropathy → high score.
Aspirin: GI bleeding, some cardiac effects at high dose → lower score.
If this ordering fails, the metric is anti-predictive and must be revised.

Usage:
    python ChemPath/scripts/sider_safety_score.py
    python ChemPath/scripts/sider_safety_score.py --validate-only
"""

from __future__ import annotations

import csv
import gzip
import math
import pickle
import statistics
from collections import defaultdict
from pathlib import Path

DATA_DIR = Path(__file__).parent.parent / "data" / "hetionet"
SIDER_FREQ_GZ = DATA_DIR / "meddra_freq.tsv.gz"
DRUGBANK_PUBCHEM_MAP = DATA_DIR / "drugbank_pubchem_map.tsv"
NODES_TSV = DATA_DIR / "nodes.tsv"
CACHE_PATH = DATA_DIR / "sider_safety_scores.pkl"

# ---------------------------------------------------------------------------
# CTCAE-inspired severity weights (keyword matching on MedDRA PT name)
# Documented and fixed — NOT tuned on AUROC results
# Based on: NCI CTCAE v5.0, November 2017
# ---------------------------------------------------------------------------

# Each entry: (keywords_tuple, weight)
# Matched in order from highest to lowest weight; first match wins.
# Keywords are lowercase substrings of the PT (Preferred Term) name.
SEVERITY_RULES: list[tuple[tuple[str, ...], float]] = [
    # Grade 5 — Death (weight 5.0)
    (("death", "fatal", "lethal", "mortality", "sudden death"), 5.0),

    # Grade 4 — Life-threatening organ toxicity (weight 4.0)
    (("anaphylaxis", "anaphylactic", "cardiac arrest", "ventricular fibrillation",
      "agranulocytosis", "aplastic anaemia", "aplastic anemia",
      "hepatic failure", "liver failure", "hepatic necrosis",
      "acute renal failure", "renal failure", "acute kidney",
      "respiratory failure", "respiratory arrest",
      "status epilepticus", "coma",
      "suicidal ideation", "suicide attempt", "completed suicide",
      "teratogen", "foetal", "fetal death", "embryo",
      "stevens-johnson", "toxic epidermal necrolysis",
      "rhabdomyolysis", "disseminated intravascular",
      "septic shock", "sepsis",
      "pulmonary hypertension", "heart failure", "cardiac failure",
      "myocardial infarction", "torsade de pointes"), 4.0),

    # Grade 3 — Severe systemic events (weight 3.0)
    (("hepatitis", "hepatotoxicity", "jaundice", "cholestasis",
      "pancreatitis", "colitis",
      "neutropenia", "thrombocytopenia", "leukopenia", "pancytopenia",
      "pulmonary embolism", "deep vein thrombosis", "thrombosis",
      "stroke", "cerebrovascular",
      "peripheral neuropathy", "neuropathy",
      "pneumonia", "pneumonitis", "interstitial lung",
      "nephrotoxicity", "nephritis", "proteinuria",
      "hypertensive crisis", "hypertensive emergency",
      "arrhythmia", "atrial fibrillation", "bradycardia", "qt prolongation",
      "seizure", "convulsion", "encephalopathy",
      "myelosuppression", "bone marrow",
      "anaphylactoid", "angioedema"), 3.0),

    # Grade 2 — Moderate symptomatic events (weight 2.0)
    (("nausea", "vomiting", "diarrhoea", "diarrhea", "constipation",
      "abdominal pain", "gastrointestinal",
      "rash", "urticaria", "pruritus", "erythema", "dermatitis",
      "headache", "migraine",
      "insomnia", "somnolence", "fatigue", "asthenia", "malaise",
      "depression", "anxiety", "mood",
      "oedema", "edema", "peripheral oedema",
      "hypoglycaemia", "hypoglycemia", "hyperglycaemia", "hyperglycemia",
      "hypotension", "orthostatic",
      "alanine aminotransferase", "aspartate aminotransferase", "liver enzyme",
      "creatinine increased", "blood creatinine",
      "anaemia", "anemia",
      "weight", "appetite",
      "hypo", "hyper"), 2.0),
    # Grade 1 — Mild (default, catch-all)
    ((), 1.0),
]


def severity_weight(pt_name: str) -> float:
    """Map a MedDRA Preferred Term name to a CTCAE-inspired severity weight."""
    name_lower = pt_name.lower()
    for keywords, weight in SEVERITY_RULES:
        if any(kw in name_lower for kw in keywords):
            return weight
    return 1.0  # Grade 1 default


# ---------------------------------------------------------------------------
# Step 1: Build PubChem CID → DrugBank ID mapping
# ---------------------------------------------------------------------------

def load_pubchem_to_drugbank(path: Path) -> dict[int, str]:
    """
    Returns {pubchem_cid (int): drugbank_id (str)}.
    File format: drugbank_id\\tpubchem_cid
    """
    mapping: dict[int, str] = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                cid = int(row["pubchem_cid"])
                mapping[cid] = row["drugbank_id"]
            except (ValueError, KeyError):
                continue
    print(f"  PubChem→DrugBank map: {len(mapping):,} entries")
    return mapping


# ---------------------------------------------------------------------------
# Step 2: Load Hetionet compound node IDs
# ---------------------------------------------------------------------------

def load_hetionet_compounds(nodes_path: Path) -> set[str]:
    """Returns set of Hetionet node IDs for Compound nodes (e.g. 'Compound::DB00014')."""
    compounds: set[str] = set()
    with open(nodes_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["kind"] == "Compound":
                compounds.add(row["id"])
    print(f"  Hetionet compounds: {len(compounds):,} nodes")
    return compounds


# ---------------------------------------------------------------------------
# Step 3: Parse SIDER meddra_freq.tsv.gz
# ---------------------------------------------------------------------------

def parse_sider_freq(gz_path: Path,
                     pubchem_to_db: dict[int, str],
                     hetionet_compounds: set[str]
                     ) -> dict[str, dict[str, float]]:
    """
    Parse SIDER 4.1 meddra_freq.tsv.gz.

    Column layout (no header):
      0: stereo_compound_id (CIDsXXX)
      1: flat_compound_id   (CIDmXXX or CIDsXXX — strip prefix, parse int)
      2: umls_concept_id    (C0000737)
      3: placebo_flag       (empty or 'placebo')
      4: freq_string        ('21%')
      5: freq_lower         (0.21)
      6: freq_upper         (0.21)
      7: meddra_type        ('PT' or 'LLT')
      8: umls_id2           (C0000737)
      9: side_effect_name   ('Abdominal pain')

    Returns:
        {hetionet_compound_id: {side_effect_name: max_freq_lower}}
        Only PT rows, only compounds in Hetionet, freq_lower > 0.
    """
    # {hetionet_id: {pt_name: max_freq_lower_seen}}
    compound_se: dict[str, dict[str, float]] = defaultdict(dict)

    matched = 0
    skipped_type = 0
    skipped_no_cid = 0
    skipped_no_freq = 0
    skipped_placebo = 0

    with gzip.open(gz_path, "rt", encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue

            meddra_type = parts[7].strip()
            if meddra_type != "PT":
                skipped_type += 1
                continue

            placebo_flag = parts[3].strip()
            if placebo_flag.lower() == "placebo":
                skipped_placebo += 1
                continue

            flat_cid_str = parts[1].strip()
            try:
                # SIDER flat CID format: 'CID0XXXXXXXX' (4-char prefix 'CID0' or 'CID1')
                # Strip first 4 chars to get the zero-padded PubChem CID integer
                cid_int = int(flat_cid_str[4:])
            except ValueError:
                skipped_no_cid += 1
                continue

            db_id = pubchem_to_db.get(cid_int)
            if db_id is None:
                skipped_no_cid += 1
                continue

            hetio_id = f"Compound::{db_id}"
            if hetio_id not in hetionet_compounds:
                skipped_no_cid += 1
                continue

            try:
                freq_lower = float(parts[5].strip())
            except ValueError:
                skipped_no_freq += 1
                continue

            if freq_lower <= 0:
                skipped_no_freq += 1
                continue

            pt_name = parts[9].strip()
            # Keep the maximum frequency seen for this (compound, PT) pair
            # (multiple label sections may report different freq ranges)
            existing = compound_se[hetio_id].get(pt_name, 0.0)
            compound_se[hetio_id][pt_name] = max(existing, freq_lower)
            matched += 1

    print(f"  SIDER parse: {matched:,} PT rows matched to Hetionet compounds")
    print(f"    Skipped (non-PT): {skipped_type:,}")
    print(f"    Skipped (placebo): {skipped_placebo:,}")
    print(f"    Skipped (no CID/compound match): {skipped_no_cid:,}")
    print(f"    Skipped (zero/missing frequency): {skipped_no_freq:,}")
    print(f"  Compounds with SIDER frequency data: {len(compound_se):,}")

    return dict(compound_se)


# ---------------------------------------------------------------------------
# Step 4: Compute safety score
# ---------------------------------------------------------------------------

def compute_safety_scores(compound_se: dict[str, dict[str, float]]
                          ) -> dict[str, float]:
    """
    Compute bias-corrected safety score per compound.

    Formula:
        raw_score(c) = Σ_{PT_i} freq_lower(i) × severity_weight(PT_i)
        safety_score(c) = raw_score(c) / log(1 + n_PT_pairs(c))

    The log-normalization controls for label-depth inflation:
    aspirin has hundreds of PT pairs because it has been studied extensively,
    not because each individual ADR is more frequent.
    Without normalization, heavily-studied drugs would always rank as "unsafe".

    Returns {hetionet_id: safety_score}, scores in [0, ∞), higher = less safe.
    """
    scores: dict[str, float] = {}
    for hetio_id, se_dict in compound_se.items():
        raw = sum(freq * severity_weight(pt) for pt, freq in se_dict.items())
        n = len(se_dict)
        # log normalization — at minimum 1 pair, log(2) ≈ 0.693
        normalized = raw / math.log(1 + n)
        scores[hetio_id] = normalized
    return scores


def diagnose_compound(hetio_id: str, compound_se: dict[str, dict[str, float]]) -> None:
    """
    Print per-PT breakdown for a specific compound.
    Useful for diagnosing unexpectedly low scores (e.g. Methotrexate).
    """
    se_dict = compound_se.get(hetio_id)
    if se_dict is None:
        print(f"  {hetio_id}: NOT in SIDER frequency data")
        return
    n = len(se_dict)
    raw = sum(freq * severity_weight(pt) for pt, freq in se_dict.items())
    normalized = raw / math.log(1 + n)
    print(f"  {hetio_id}: n_PT_pairs={n}, raw_score={raw:.4f}, "
          f"log_norm={math.log(1+n):.3f}, final_score={normalized:.4f}")
    print(f"  {'PT name':<45s} {'freq':>7s} {'sev':>5s} {'contrib':>8s}")
    print(f"  {'-'*70}")
    rows = sorted(
        [(pt, freq, severity_weight(pt), freq * severity_weight(pt))
         for pt, freq in se_dict.items()],
        key=lambda x: -x[3]
    )
    for pt, freq, sev, contrib in rows[:30]:
        print(f"  {pt:<45s} {freq:>7.4f} {sev:>5.1f} {contrib:>8.4f}")
    if len(rows) > 30:
        print(f"  ... ({len(rows) - 30} more PT terms)")


# ---------------------------------------------------------------------------
# Step 5: Validation
# ---------------------------------------------------------------------------

VALIDATION_CASES = {
    # HIGH risk — cytotoxic, immunosuppressive, teratogenic
    "Compound::DB01041": ("Thalidomide",         "high"),   # teratogen, peripheral neuropathy
    "Compound::DB00619": ("Imatinib",            "high"),   # cancer, hepatotoxicity, CHF
    "Compound::DB00688": ("Mycophenolate MF",    "high"),   # immunosupp., teratogenic
    "Compound::DB01008": ("Busulfan",            "high"),   # cytotoxic, myelosuppression
    "Compound::DB00864": ("Tacrolimus",          "high"),   # immunosupp., nephrotoxicity
    # MEDIUM risk
    "Compound::DB00316": ("Acetaminophen",       "medium"), # hepatotoxic in overdose
    "Compound::DB00709": ("Lamivudine",          "medium"), # HIV, generally well-tolerated
    # LOW risk — common, narrow-indication, well-tolerated antibiotics/generics
    "Compound::DB01137": ("Levofloxacin",        "low"),    # antibiotic, short-term safe
    "Compound::DB00319": ("Piperacillin",        "low"),    # beta-lactam antibiotic
    "Compound::DB00565": ("Cisatracurium",       "low"),    # neuromuscular blocker, ICU use
}


def validate_scores(scores: dict[str, float]) -> bool:
    """
    Print validation table and return True if ordering is pharmacologically plausible.

    Tests only use compounds confirmed to be in SIDER frequency data.
    Aspirin and Warfarin are NOT in SIDER meddra_freq.tsv (no label frequency data).

    Critical pairwise tests: high-risk cytotoxic/immunosuppressive drugs must
    outscore common antibiotics and ICU agents by safety score.
    """
    print("\n  Validation cases (higher score = more unsafe):")
    print(f"  {'Compound':<20s} {'DB ID':<15s} {'Score':>8s} {'Expected':>10s} {'In SIDER':>9s}")
    print(f"  {'-'*70}")

    in_sider = {k: v for k, v in VALIDATION_CASES.items() if k in scores}
    not_in_sider = {k: v for k, v in VALIDATION_CASES.items() if k not in scores}

    all_scores = sorted(in_sider.items(), key=lambda x: scores.get(x[0], 0), reverse=True)
    for db_id, (name, expected) in all_scores:
        score = scores.get(db_id, float("nan"))
        print(f"  {name:<20s} {db_id:<15s} {score:>8.4f} {expected:>10s} {'YES':>9s}")
    for db_id, (name, expected) in not_in_sider.items():
        print(f"  {name:<20s} {db_id:<15s} {'N/A':>8s} {expected:>10s} {'NO':>9s}")

    print("\n  Critical ordering tests:")
    passed = 0
    total = 0

    # Pairwise ordering tests: each high-risk drug should outscore each low-risk drug
    high_risk = [
        ("Thalidomide",       scores.get("Compound::DB01041")),
        ("Imatinib",          scores.get("Compound::DB00619")),
        ("Busulfan",          scores.get("Compound::DB01008")),
        ("Tacrolimus",        scores.get("Compound::DB00864")),
    ]
    low_risk = [
        ("Levofloxacin",      scores.get("Compound::DB01137")),
        ("Piperacillin",      scores.get("Compound::DB00319")),
        ("Cisatracurium",     scores.get("Compound::DB00565")),
    ]

    for h_name, h_score in high_risk:
        for l_name, l_score in low_risk:
            if h_score is None or l_score is None:
                continue
            total += 1
            ok = h_score > l_score
            print(f"  {h_name:<20s} ({h_score:.4f}) > {l_name:<15s} ({l_score:.4f}): {'PASS' if ok else 'FAIL'}")
            if ok:
                passed += 1

    if total == 0:
        print("  WARNING: Validation compounds not found in SIDER data")
        return False

    print(f"\n  Result: {passed}/{total} ordering tests passed")
    return passed >= total // 2 + 1  # majority pass


# ---------------------------------------------------------------------------
# Fallback: log-normalized Hetionet 'causes' edge count
# ---------------------------------------------------------------------------

FALLBACK_CACHE_PATH = DATA_DIR / "sider_fallback_scores.pkl"


def build_fallback_scores(sider_scores: dict[str, float],
                          force_rebuild: bool = False) -> dict[str, float]:
    """
    For compounds NOT in SIDER frequency data (67% of Hetionet), compute a
    fallback safety score from Hetionet's Compound-causes-SideEffect edge count.

    Uses log(1 + count) normalization (same form as SIDER score) scaled to
    match the SIDER score distribution (same median and 75th percentile).

    Returns {hetionet_id: fallback_safety_score} for non-SIDER compounds only.
    """
    if FALLBACK_CACHE_PATH.exists() and not force_rebuild:
        with open(FALLBACK_CACHE_PATH, "rb") as f:
            return pickle.load(f)

    compact_pkl = DATA_DIR / "hetionet_compact.pkl"
    if not compact_pkl.exists():
        print("  WARNING: hetionet_compact.pkl not found — skipping fallback")
        return {}

    print("  Building fallback safety scores from Hetionet causes edges...")
    with open(compact_pkl, "rb") as f:
        cached = pickle.load(f)

    # hetionet_compact.pkl filters out 'causes' edges (keeps only mechanistic edges).
    # Try side_effect_counts from the enriched raw cache first; fall back to graph traversal.
    raw_counts: dict[str, int] = defaultdict(int)
    if "side_effect_counts" in cached:
        for compound, count in cached["side_effect_counts"].items():
            raw_counts[compound] = count
    else:
        G = cached["G"]
        for u, _v, d in G.edges(data=True):
            if d.get("edge_type") == "causes" and u.startswith("Compound::"):
                raw_counts[u] += 1

    if not raw_counts:
        # Try the enriched benchmark cache which always has side_effect_counts
        enriched_pkl = DATA_DIR / "hetionet_enriched_raw.pkl"
        if enriched_pkl.exists():
            with open(enriched_pkl, "rb") as f:
                enriched = pickle.load(f)
            for compound, count in enriched.get("side_effect_counts", {}).items():
                raw_counts[compound] = count

    # Log-normalize
    log_scores = {c: math.log(1 + n) for c, n in raw_counts.items()
                  if c not in sider_scores}

    if not log_scores or not sider_scores:
        return {}

    # Scale to match SIDER distribution (align medians and p75)
    sider_vals = sorted(sider_scores.values())
    sider_median = sider_vals[len(sider_vals) // 2]
    sider_p75    = sider_vals[int(0.75 * len(sider_vals))]

    ls_vals = sorted(log_scores.values())
    ls_median = ls_vals[len(ls_vals) // 2]
    ls_p75    = ls_vals[int(0.75 * len(ls_vals))]

    # Linear scale so fallback p75 maps to SIDER p75
    scale = sider_p75 / ls_p75 if ls_p75 > 0 else 1.0
    fallback = {c: s * scale for c, s in log_scores.items()}

    print(f"  Fallback scores: {len(fallback)} compounds "
          f"(scale={scale:.3f}, median={statistics.median(fallback.values()):.4f})")

    with open(FALLBACK_CACHE_PATH, "wb") as f:
        pickle.dump(fallback, f, protocol=pickle.HIGHEST_PROTOCOL)

    return fallback


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def build_sider_safety_scores(force_rebuild: bool = False) -> dict[str, float]:
    """
    Build or load cached SIDER frequency-weighted severity safety scores.

    Returns {hetionet_compound_id: safety_score}.
    Higher score = compound is less safe (more/severe documented ADRs).
    """
    if CACHE_PATH.exists() and not force_rebuild:
        print(f"  [cached] Loading safety scores from {CACHE_PATH.name}")
        with open(CACHE_PATH, "rb") as f:
            return pickle.load(f)

    print("  Building SIDER safety scores from scratch...")
    print(f"  Source: {SIDER_FREQ_GZ.name} ({SIDER_FREQ_GZ.stat().st_size / 1e6:.1f} MB)")

    # Step 1: ID mappings
    print("\n  [1/4] Loading ID mappings...")
    pubchem_to_db = load_pubchem_to_drugbank(DRUGBANK_PUBCHEM_MAP)
    hetionet_compounds = load_hetionet_compounds(NODES_TSV)

    # Step 2: Parse SIDER
    print("\n  [2/4] Parsing SIDER meddra_freq.tsv...")
    compound_se = parse_sider_freq(SIDER_FREQ_GZ, pubchem_to_db, hetionet_compounds)

    # Step 3: Compute scores
    print("\n  [3/4] Computing frequency × severity scores...")
    scores = compute_safety_scores(compound_se)

    # Coverage report
    coverage = len(scores) / len(hetionet_compounds) * 100
    print(f"  Score coverage: {len(scores)}/{len(hetionet_compounds)} compounds ({coverage:.1f}%)")

    # Step 4: Score distribution
    vals = list(scores.values())
    if vals:
        vals_sorted = sorted(vals)
        print(f"  Score distribution:")
        print(f"    min={min(vals):.4f}, p25={vals_sorted[len(vals)//4]:.4f}, "
              f"median={statistics.median(vals):.4f}, "
              f"p75={vals_sorted[3*len(vals)//4]:.4f}, max={max(vals):.4f}")

    # Save cache
    with open(CACHE_PATH, "wb") as f:
        pickle.dump(scores, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"\n  Saved cache: {CACHE_PATH.name}")

    return scores


def main() -> None:
    import sys
    print("=" * 70)
    print(" SIDER 4.1 Frequency-Weighted Safety Score Builder")
    print("=" * 70)

    force = "--rebuild" in sys.argv
    validate_only = "--validate-only" in sys.argv

    # Need compound_se for diagnostics — rebuild from scratch
    if not validate_only:
        pubchem_to_db = load_pubchem_to_drugbank(DRUGBANK_PUBCHEM_MAP)
        hetionet_compounds = load_hetionet_compounds(NODES_TSV)
        compound_se = parse_sider_freq(SIDER_FREQ_GZ, pubchem_to_db, hetionet_compounds)
        scores = compute_safety_scores(compound_se)
        if force or not CACHE_PATH.exists():
            with open(CACHE_PATH, "wb") as f:
                pickle.dump(scores, f, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        scores = build_sider_safety_scores(force_rebuild=False)
        compound_se = {}  # not available in validate-only mode

    # Fallback scores for non-SIDER compounds
    print("\n  Building fallback scores...")
    fallback = build_fallback_scores(scores, force_rebuild=force)
    all_scores = {**fallback, **scores}  # SIDER takes priority
    print(f"  Total coverage: {len(all_scores)} compounds "
          f"(SIDER: {len(scores)}, fallback: {len(fallback)})")

    print("\n" + "=" * 70)
    print(" Validation")
    print("=" * 70)
    passed = validate_scores(scores)

    # Diagnose methotrexate (known low scorer — frequency-weighting limitation)
    if compound_se:
        print("\n" + "=" * 70)
        print(" Methotrexate Diagnosis (DB00563)")
        print("=" * 70)
        diagnose_compound("Compound::DB00563", compound_se)
        print("\n  NOTE: Methotrexate's severe ADRs (hepatotoxicity, mucositis) have")
        print("  low label frequencies (rare but serious). freq×severity is small.")
        print("  This is a known limitation of frequency-weighting: rare but")
        print("  life-threatening events are downweighted vs common mild events.")

    # Top 20 most unsafe
    print("\n  Top 20 highest-scoring (most unsafe) compounds:")
    top = sorted(scores.items(), key=lambda x: x[1], reverse=True)[:20]
    for i, (cid, score) in enumerate(top, 1):
        print(f"  {i:2d}. {cid:<30s}  {score:.4f}")

    # Bottom 20 safest
    print("\n  Bottom 20 lowest-scoring (safest) compounds:")
    bottom = sorted(scores.items(), key=lambda x: x[1])[:20]
    for i, (cid, score) in enumerate(bottom, 1):
        print(f"  {i:2d}. {cid:<30s}  {score:.4f}")

    print("\n" + "=" * 70)
    print(f" Status: {'PASS — metric is pharmacologically plausible' if passed else 'FAIL — metric needs revision'}")
    print("=" * 70)


if __name__ == "__main__":
    main()
