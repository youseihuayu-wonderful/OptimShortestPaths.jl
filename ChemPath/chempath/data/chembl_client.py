"""
ChEMBL REST API client — direct HTTP, no external dependencies.

Fetches bioactivity (IC50) data, compound metadata, and target info
from ChEMBL's public REST API (https://www.ebi.ac.uk/chembl/api/data/).

Features:
  - Paginated fetching with configurable limits
  - Local JSON cache to avoid redundant API calls
  - Rate limiting (polite 0.5s delay between requests)
"""

import json
import time
import urllib.request
import urllib.error
from pathlib import Path

BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"
CACHE_DIR = Path(__file__).parent.parent.parent / ".chembl_cache"
REQUEST_DELAY = 0.5  # seconds between API calls


def _fetch_json(url: str, use_cache: bool = True) -> dict:
    """Fetch JSON from URL with optional file-based caching."""
    if use_cache:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        cache_key = url.replace("/", "_").replace(":", "_").replace("?", "_").replace("&", "_")
        # Truncate long cache keys
        if len(cache_key) > 200:
            import hashlib
            cache_key = hashlib.md5(url.encode()).hexdigest()
        cache_file = CACHE_DIR / f"{cache_key}.json"
        if cache_file.exists():
            return json.loads(cache_file.read_text())

    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode())
    except urllib.error.HTTPError as e:
        raise RuntimeError(f"ChEMBL API error {e.code}: {e.reason} for {url}") from e
    except urllib.error.URLError as e:
        raise RuntimeError(f"Network error: {e.reason} for {url}") from e

    if use_cache:
        cache_file.write_text(json.dumps(data))

    time.sleep(REQUEST_DELAY)
    return data


def _fetch_all_pages(url: str, key: str, max_records: int = 1000, use_cache: bool = True) -> list:
    """Fetch all pages of a paginated ChEMBL endpoint."""
    results = []
    offset = 0
    limit = min(100, max_records)

    while len(results) < max_records:
        page_url = f"{url}&offset={offset}&limit={limit}"
        data = _fetch_json(page_url, use_cache=use_cache)
        page_results = data.get(key, [])
        if not page_results:
            break
        results.extend(page_results)
        offset += limit

        meta = data.get("page_meta", {})
        if not meta.get("next"):
            break

    return results[:max_records]


def fetch_activities_for_target(
    target_chembl_id: str,
    standard_type: str = "IC50",
    max_records: int = 500,
    use_cache: bool = True,
) -> list[dict]:
    """
    Fetch bioactivity records for a target from ChEMBL.

    Returns list of dicts with keys:
      molecule_chembl_id, molecule_pref_name, canonical_smiles,
      standard_value, standard_units, standard_type,
      target_chembl_id, target_pref_name, pchembl_value,
      assay_type, data_validity_comment
    """
    url = (
        f"{BASE_URL}/activity.json?"
        f"target_chembl_id={target_chembl_id}"
        f"&standard_type={standard_type}"
        f"&standard_relation=%3D"  # only exact "=" measurements
    )
    return _fetch_all_pages(url, "activities", max_records=max_records, use_cache=use_cache)


def fetch_molecule(chembl_id: str, use_cache: bool = True) -> dict:
    """Fetch molecule metadata (max_phase, mol weight, etc.)."""
    url = f"{BASE_URL}/molecule/{chembl_id}.json"
    return _fetch_json(url, use_cache=use_cache)


def fetch_target(chembl_id: str, use_cache: bool = True) -> dict:
    """Fetch target metadata."""
    url = f"{BASE_URL}/target/{chembl_id}.json"
    return _fetch_json(url, use_cache=use_cache)


# --- Anti-cancer targets for Phase 2 ---

ANTICANCER_TARGETS = [
    {"chembl_id": "CHEMBL203", "name": "EGFR", "description": "Epidermal Growth Factor Receptor"},
    {"chembl_id": "CHEMBL1824", "name": "HER2", "description": "Receptor tyrosine-protein kinase erbB-2"},
    {"chembl_id": "CHEMBL4282", "name": "ALK", "description": "ALK tyrosine kinase receptor"},
    {"chembl_id": "CHEMBL5145", "name": "BRAF", "description": "Serine/threonine-protein kinase B-raf"},
    {"chembl_id": "CHEMBL1862", "name": "ABL1", "description": "Tyrosine-protein kinase ABL1"},
    {"chembl_id": "CHEMBL2842", "name": "VEGFR2", "description": "Vascular endothelial growth factor receptor 2"},
    {"chembl_id": "CHEMBL4005", "name": "MET", "description": "Hepatocyte growth factor receptor"},
    {"chembl_id": "CHEMBL3594", "name": "FGFR1", "description": "Fibroblast growth factor receptor 1"},
    {"chembl_id": "CHEMBL4630", "name": "PI3Kα", "description": "PI3-kinase p110-alpha subunit"},
    {"chembl_id": "CHEMBL2185", "name": "mTOR", "description": "Serine/threonine-protein kinase mTOR"},
]


def clean_activity(raw: dict) -> dict | None:
    """
    Clean a single raw ChEMBL activity record into our standard format.
    Returns None if the record is unusable.
    """
    value = raw.get("standard_value")
    units = raw.get("standard_units")
    smiles = raw.get("canonical_smiles")
    compound_id = raw.get("molecule_chembl_id")
    target_id = raw.get("target_chembl_id")

    # Skip records with missing critical fields
    if value is None or smiles is None or compound_id is None or target_id is None:
        return None

    try:
        value = float(value)
    except (ValueError, TypeError):
        return None

    if value <= 0:
        return None

    # Skip flagged data
    if raw.get("data_validity_comment"):
        return None

    # Standardize units to nM
    if units == "nM":
        ic50_nm = value
    elif units == "uM":
        ic50_nm = value * 1000.0
    elif units == "pM":
        ic50_nm = value / 1000.0
    elif units == "M":
        ic50_nm = value * 1e9
    else:
        return None  # unknown unit

    return {
        "compound": compound_id,
        "compound_name": raw.get("molecule_pref_name") or compound_id,
        "smiles": smiles,
        "target": target_id,
        "target_name": raw.get("target_pref_name") or target_id,
        "type": raw.get("standard_type", "IC50"),
        "value": round(ic50_nm, 2),
        "units": "nM",
        "source": "experimental",
        "pchembl_value": raw.get("pchembl_value"),
        "assay_type": raw.get("assay_type"),
    }


def fetch_and_clean_target_data(
    target_chembl_id: str,
    max_records: int = 500,
    use_cache: bool = True,
    verbose: bool = True,
) -> list[dict]:
    """Fetch and clean IC50 data for a single target."""
    if verbose:
        print(f"  Fetching IC50 data for {target_chembl_id}...")

    raw = fetch_activities_for_target(
        target_chembl_id, max_records=max_records, use_cache=use_cache
    )
    if verbose:
        print(f"    Raw records: {len(raw)}")

    cleaned = []
    for r in raw:
        c = clean_activity(r)
        if c is not None:
            cleaned.append(c)

    # Deduplicate: keep the record with the lowest IC50 per compound-target pair
    best = {}
    for c in cleaned:
        key = (c["compound"], c["target"])
        if key not in best or c["value"] < best[key]["value"]:
            best[key] = c
    deduped = list(best.values())

    if verbose:
        print(f"    After cleaning: {len(cleaned)}, after dedup: {len(deduped)}")

    return deduped


def fetch_real_data(
    targets: list[dict] | None = None,
    max_per_target: int = 500,
    use_cache: bool = True,
    verbose: bool = True,
) -> dict:
    """
    Fetch real ChEMBL data for multiple targets.
    Returns data in the same format as load_mock_data().
    """
    if targets is None:
        targets = ANTICANCER_TARGETS

    if verbose:
        print(f"Fetching data for {len(targets)} targets from ChEMBL...")

    all_activities = []
    compounds_seen = {}
    targets_out = []

    for t in targets:
        tid = t["chembl_id"]
        targets_out.append({
            "chembl_id": tid,
            "name": t["name"],
            "organism": "Homo sapiens",
            "type": "SINGLE PROTEIN",
        })

        activities = fetch_and_clean_target_data(
            tid, max_records=max_per_target, use_cache=use_cache, verbose=verbose
        )

        for a in activities:
            cid = a["compound"]
            if cid not in compounds_seen:
                compounds_seen[cid] = {
                    "chembl_id": cid,
                    "name": a["compound_name"],
                    "smiles": a["smiles"],
                    "phase": 0,  # will be enriched later if needed
                }
            all_activities.append({
                "compound": cid,
                "target": a["target"],
                "type": "IC50",
                "value": a["value"],
                "units": "nM",
                "source": "experimental",
            })

    compounds_out = list(compounds_seen.values())

    if verbose:
        print(f"\nSummary: {len(compounds_out)} unique compounds, "
              f"{len(targets_out)} targets, {len(all_activities)} bioactivities")

    return {
        "compounds": compounds_out,
        "targets": targets_out,
        "bioactivities": all_activities,
        "toxicity": {},  # no toxicity data from ChEMBL directly
    }


def enrich_compounds_with_phase(
    data: dict,
    batch_size: int = 50,
    use_cache: bool = True,
    verbose: bool = True,
) -> dict:
    """
    Enrich compound records with max_phase from ChEMBL molecule endpoint.
    Modifies data in place and returns it.
    """
    compounds = data["compounds"]
    total = len(compounds)
    enriched = 0

    if verbose:
        print(f"Enriching {total} compounds with clinical phase data...")

    for i in range(0, total, batch_size):
        batch = compounds[i:i + batch_size]
        ids = ",".join(c["chembl_id"] for c in batch)
        url = f"{BASE_URL}/molecule.json?molecule_chembl_id__in={ids}&limit={batch_size}"

        try:
            resp = _fetch_json(url, use_cache=use_cache)
            molecules = resp.get("molecules", [])
            phase_map = {}
            for mol in molecules:
                mid = mol.get("molecule_chembl_id")
                phase = mol.get("max_phase")
                if mid and phase is not None:
                    try:
                        phase_map[mid] = int(float(phase))
                    except (ValueError, TypeError):
                        phase_map[mid] = 0

            for c in batch:
                if c["chembl_id"] in phase_map:
                    c["phase"] = phase_map[c["chembl_id"]]
                    if phase_map[c["chembl_id"]] > 0:
                        enriched += 1
        except RuntimeError as e:
            if verbose:
                print(f"    Warning: batch {i//batch_size + 1} failed: {e}")
            continue

        if verbose and (i + batch_size) % 500 == 0:
            print(f"    Processed {min(i + batch_size, total)}/{total}...")

    if verbose:
        print(f"  Enriched {enriched}/{total} compounds with clinical phase > 0")

    return data


def save_data(data: dict, path: str | Path) -> None:
    """Save fetched data to JSON for offline use."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2))
    print(f"Saved to {path} ({path.stat().st_size / 1024:.1f} KB)")


def load_saved_data(path: str | Path) -> dict:
    """Load previously saved data from JSON."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"No saved data at {path}")
    return json.loads(path.read_text())
