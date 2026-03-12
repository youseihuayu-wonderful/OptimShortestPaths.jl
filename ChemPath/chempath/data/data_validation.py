"""
Automated data validation for curated dataset.

Catches mislabeling, duplicate IDs, orphaned references, and out-of-range
values before they propagate into the graph or benchmark.

Usage:
    from chempath.data.data_validation import validate_dataset
    errors = validate_dataset()          # returns list of error strings
    assert not errors, "\\n".join(errors)  # fail fast if invalid

Called automatically by load_curated_data() in debug mode.
"""

from __future__ import annotations


# Known ChEMBL ID → name mapping for compounds in our dataset.
# Source: ChEMBL 36 database.  Used to catch mislabeling (e.g. CHEMBL1336
# was once incorrectly labeled "Lapatinib" instead of "Sorafenib").
CHEMBL_COMPOUND_NAMES = {
    "CHEMBL939": "Gefitinib",
    "CHEMBL553": "Erlotinib",
    "CHEMBL1201585": "Osimertinib",
    "CHEMBL554": "Lapatinib",
    "CHEMBL180022": "Neratinib",
    "CHEMBL24828": "Vandetanib",
    "CHEMBL3545252": "Alectinib",
    "CHEMBL2007641": "Ceritinib",
    "CHEMBL601719": "Crizotinib",
    "CHEMBL1229517": "Vemurafenib",
    "CHEMBL2028663": "Dabrafenib",
    "CHEMBL941": "Imatinib",
    "CHEMBL288441": "Bosutinib",
    "CHEMBL255863": "Nilotinib",
    "CHEMBL1171837": "Ponatinib",
    "CHEMBL1421": "Dasatinib",
    "CHEMBL1336": "Sorafenib",
    "CHEMBL1642": "Everolimus",
}

# Known ChEMBL ID → name mapping for targets.
CHEMBL_TARGET_NAMES = {
    "CHEMBL203": "EGFR",
    "CHEMBL1824": "HER2",
    "CHEMBL4282": "ALK",
    "CHEMBL5145": "BRAF",
    "CHEMBL1862": "ABL1",
    "CHEMBL2842": "VEGFR2",
    "CHEMBL4005": "MET",
    "CHEMBL3594": "FGFR1",
    "CHEMBL4630": "PI3Ka",
    "CHEMBL2185": "mTOR",
    "CHEMBL1936": "KIT",
}

# Plausible IC50 range (nM).  Values outside this trigger a warning.
IC50_MIN = 0.01    # sub-picomolar — suspiciously potent
IC50_MAX = 100_000  # 100 µM — essentially inactive


def validate_dataset(data: dict | None = None) -> list[str]:
    """
    Run all validation checks on the curated dataset.

    Parameters
    ----------
    data : dict or None
        Dataset dict with keys: compounds, targets, bioactivities, toxicity.
        If None, loads from curated_data.py automatically.

    Returns
    -------
    list[str]
        List of error messages.  Empty list means all checks passed.
    """
    if data is None:
        from chempath.data.curated_data import load_curated_data
        data = load_curated_data()

    errors: list[str] = []
    errors.extend(_check_unique_ids(data))
    errors.extend(_check_compound_names(data))
    errors.extend(_check_target_names(data))
    errors.extend(_check_bioactivity_references(data))
    errors.extend(_check_ic50_range(data))
    errors.extend(_check_duplicate_edges(data))
    errors.extend(_check_toxicity_coverage(data))
    errors.extend(_check_ground_truth(data))
    return errors


def _check_unique_ids(data: dict) -> list[str]:
    """Compound and target ChEMBL IDs must be unique."""
    errors = []
    cids = [c["chembl_id"] for c in data["compounds"]]
    seen = set()
    for cid in cids:
        if cid in seen:
            errors.append(f"Duplicate compound ID: {cid}")
        seen.add(cid)

    tids = [t["chembl_id"] for t in data["targets"]]
    seen = set()
    for tid in tids:
        if tid in seen:
            errors.append(f"Duplicate target ID: {tid}")
        seen.add(tid)
    return errors


def _check_compound_names(data: dict) -> list[str]:
    """Verify compound names match known ChEMBL mappings."""
    errors = []
    for c in data["compounds"]:
        cid = c["chembl_id"]
        if cid in CHEMBL_COMPOUND_NAMES:
            expected = CHEMBL_COMPOUND_NAMES[cid]
            if c["name"] != expected:
                errors.append(
                    f"Compound mislabel: {cid} is '{c['name']}' "
                    f"but should be '{expected}'"
                )
    return errors


def _check_target_names(data: dict) -> list[str]:
    """Verify target names match known ChEMBL mappings."""
    errors = []
    for t in data["targets"]:
        tid = t["chembl_id"]
        if tid in CHEMBL_TARGET_NAMES:
            expected = CHEMBL_TARGET_NAMES[tid]
            if t["name"] != expected:
                errors.append(
                    f"Target mislabel: {tid} is '{t['name']}' "
                    f"but should be '{expected}'"
                )
    return errors


def _check_bioactivity_references(data: dict) -> list[str]:
    """Every bioactivity must reference a valid compound and target."""
    errors = []
    cid_set = {c["chembl_id"] for c in data["compounds"]}
    tid_set = {t["chembl_id"] for t in data["targets"]}

    for i, b in enumerate(data["bioactivities"]):
        if b["compound"] not in cid_set:
            errors.append(
                f"Bioactivity[{i}]: unknown compound '{b['compound']}'"
            )
        if b["target"] not in tid_set:
            errors.append(
                f"Bioactivity[{i}]: unknown target '{b['target']}'"
            )
    return errors


def _check_ic50_range(data: dict) -> list[str]:
    """IC50 values should be within a plausible range."""
    errors = []
    for i, b in enumerate(data["bioactivities"]):
        val = b.get("value")
        if val is None:
            errors.append(f"Bioactivity[{i}]: missing IC50 value")
            continue
        if val <= 0:
            errors.append(
                f"Bioactivity[{i}]: IC50={val} nM is non-positive "
                f"({b['compound']}→{b['target']})"
            )
        elif val < IC50_MIN:
            errors.append(
                f"Bioactivity[{i}]: IC50={val} nM is suspiciously low "
                f"({b['compound']}→{b['target']})"
            )
        elif val > IC50_MAX:
            errors.append(
                f"Bioactivity[{i}]: IC50={val} nM is above {IC50_MAX} "
                f"({b['compound']}→{b['target']})"
            )
    return errors


def _check_duplicate_edges(data: dict) -> list[str]:
    """Each (compound, target) pair should appear at most once."""
    errors = []
    seen = set()
    for i, b in enumerate(data["bioactivities"]):
        pair = (b["compound"], b["target"])
        if pair in seen:
            errors.append(
                f"Bioactivity[{i}]: duplicate edge {pair[0]}→{pair[1]}"
            )
        seen.add(pair)
    return errors


def _check_toxicity_coverage(data: dict) -> list[str]:
    """Every non-test compound should have toxicity data, and vice versa."""
    errors = []
    valid_cids = {
        c["chembl_id"]
        for c in data["compounds"]
        if not c["chembl_id"].startswith("CHEMBL_")
    }
    tox_ids = set(data.get("toxicity", {}).keys())

    missing = valid_cids - tox_ids
    for cid in sorted(missing):
        errors.append(f"Compound {cid} has no toxicity data")

    extra = tox_ids - valid_cids
    for cid in sorted(extra):
        errors.append(f"Toxicity data for unknown compound {cid}")
    return errors


def _check_ground_truth(data: dict) -> list[str]:
    """Ground truth compounds must exist in the dataset."""
    errors = []
    try:
        from chempath.graph.benchmark import GROUND_TRUTH
    except ImportError:
        return errors

    cid_set = {c["chembl_id"] for c in data["compounds"]}
    for cid in GROUND_TRUTH:
        if cid not in cid_set:
            errors.append(
                f"Ground truth compound {cid} not found in dataset"
            )
    return errors
