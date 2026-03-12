"""
ChemPath POC v3: Multi-Objective Pareto Drug Repurposing Benchmark.

Key fixes over POC2:
  - Efficacy: pIC50 modulates ±30% around base edge-type weight (preserves topology)
  - Safety: SIDER 4.1 freq×CTCAE-severity with log-normalization (sider_safety_score.py)
  - Evidence: unchanged from POC2

Enriches Hetionet's binary graph with:
  - Continuous efficacy weights from ChEMBL pIC50 bioactivity data
  - Severity-weighted safety scores from SIDER frequency-tier data
  - Per-edge evidence scores from Hetionet internal metadata

Then benchmarks multi-objective scoring against PharmacotherapyDB.

Usage:
    python chempath_enriched_benchmark.py [--skip-chembl] [--sensitivity]
"""
from __future__ import annotations

import bz2
import csv
import gzip
import io
import json
import math
import os
import pickle
import random
import sys
import time
import urllib.error
import urllib.request
from collections import defaultdict
from pathlib import Path

import networkx as nx

# Kim's SIDER 4.1 bias-corrected safety score (freq × CTCAE severity / log norm)
sys.path.insert(0, str(Path(__file__).parent))
from sider_safety_score import build_sider_safety_scores, build_fallback_scores

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
HETIONET_JSON_URL = (
    "https://github.com/hetio/hetionet/raw/main/hetnet/json/hetionet-v1.0.json.bz2"
)
PHARMACOTHERAPYDB_URL = (
    "https://raw.githubusercontent.com/dhimmel/indications/"
    "gh-pages/catalog/indications.tsv"
)
DRUGBANK_CHEMBL_MAP_URL = (
    "https://raw.githubusercontent.com/dhimmel/drugbank/"
    "gh-pages/data/mapping/chembl.tsv"
)
DRUGBANK_PUBCHEM_MAP_URL = (
    "https://raw.githubusercontent.com/dhimmel/drugbank/"
    "3e87872db5fca5ac427ce27464ab945c0ceb4ec6/data/pubchem-mapping.tsv"
)
SIDER_FREQ_URL = (
    "https://raw.githubusercontent.com/dhimmel/SIDER4/"
    "master/download/meddra_freq.tsv.gz"
)
CHEMBL_ACTIVITY_URL = (
    "https://www.ebi.ac.uk/chembl/api/data/activity.json"
    "?molecule_chembl_id={mol}&pchembl_value__isnull=false"
    "&limit=1000&format=json"
)

CACHE_DIR = Path(__file__).resolve().parent.parent / "data" / "hetionet"
CHEMBL_CACHE_DIR = CACHE_DIR / "chembl_cache"

RELEVANT_EDGES = {
    "binds", "downregulates", "upregulates", "interacts",
    "associates", "palliates", "resembles", "includes",
}

# Edge-type base probabilities (used when no enrichment data available)
EDGE_TYPE_PROB = {
    "binds": 0.80, "downregulates": 0.75, "upregulates": 0.75,
    "interacts": 0.70, "associates": 0.80, "palliates": 0.60,
    "resembles": 0.50, "includes": 0.60, "regulates": 0.50,
    "expresses": 0.50, "localizes": 0.70,
}

# SIDER frequency tier -> severity weight (higher = more clinically concerning)
FREQ_TIER_WEIGHTS = {
    "postmarketing": 0.1,
    "rare": 0.3,
    "infrequent": 0.5,
    "frequent": 0.8,
    "common": 0.8,
    "very common": 1.0,
    "very frequent": 1.0,
}

# Modulation range for pIC50 adjustment (±MODULATION_RANGE around base weight)
# Sensitivity: tested 0.30 (default) and 0.50 — 0.30 is optimal
MODULATION_RANGE = 0.30


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------
def _download(url: str, dest: Path, desc: str = "") -> Path:
    if dest.exists():
        print(f"  [cached] {dest.name}")
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"  Downloading {desc or dest.name}...")
    try:
        urllib.request.urlretrieve(url, dest)
    except Exception as e:
        print(f"  WARN: download failed ({e}), trying fallback...")
        raise
    return dest


def _sigmoid(x: float, center: float = 6.0) -> float:
    """Sigmoid centered at `center`, maps (-inf, inf) -> (0, 1)."""
    z = x - center
    if z > 20:
        return 1.0
    if z < -20:
        return 0.0
    return 1.0 / (1.0 + math.exp(-z))


def _safe_log(p: float) -> float:
    """Convert probability to additive cost: -log(p), clamped."""
    p = max(p, 1e-6)
    p = min(p, 1.0 - 1e-6)
    return -math.log(p)


# ---------------------------------------------------------------------------
# Part 0: Load Hetionet with FULL edge metadata
# ---------------------------------------------------------------------------
def load_hetionet_enriched() -> dict:
    """Load Hetionet with edge metadata preserved for enrichment."""
    enriched_pkl = CACHE_DIR / "hetionet_enriched_raw.pkl"
    if enriched_pkl.exists():
        print("  [cached] enriched raw data")
        with open(enriched_pkl, "rb") as f:
            return pickle.load(f)

    # Download raw JSON
    bz2_path = _download(HETIONET_JSON_URL, CACHE_DIR / "hetionet-v1.0.json.bz2")
    print("  Decompressing and parsing JSON...")
    with bz2.open(bz2_path, "rt") as f:
        raw = json.load(f)

    # Build node lookup
    node_map = {}  # (kind, identifier) -> hetionet_id
    for n in raw["nodes"]:
        het_id = f"{n['kind']}::{n['identifier']}"
        node_map[(n["kind"], str(n["identifier"]))] = het_id

    # Build graph with metadata
    G = nx.DiGraph()
    for n in raw["nodes"]:
        het_id = f"{n['kind']}::{n['identifier']}"
        G.add_node(het_id, kind=n["kind"], name=n.get("name", ""))

    ground_truth = defaultdict(list)  # compound -> [diseases]
    side_effect_counts = defaultdict(int)
    edge_metadata = {}  # (u, v) -> metadata dict

    for e in raw["edges"]:
        kind = e["kind"]
        src_key = (e["source_id"][0], str(e["source_id"][1]))
        tgt_key = (e["target_id"][0], str(e["target_id"][1]))
        src_id = node_map.get(src_key)
        tgt_id = node_map.get(tgt_key)
        if not src_id or not tgt_id:
            continue

        data = e.get("data", {})

        # Hold out treats edges
        if kind == "treats":
            ground_truth[src_id].append(tgt_id)
            continue

        # Count side effects
        if kind == "causes" and src_id.startswith("Compound::"):
            side_effect_counts[src_id] += 1
            continue

        # Filter to relevant edges
        if kind not in RELEVANT_EDGES:
            continue

        # Store metadata for enrichment
        meta = {
            "edge_type": kind,
            "pubmed_ids": data.get("pubmed_ids", []),
            "sources": data.get("sources", []) if isinstance(data.get("sources"), list) else ([data["source"]] if "source" in data else []),
            "actions": data.get("actions", []),
            "similarity": data.get("similarity"),
            "unbiased": data.get("unbiased", False),
        }

        direction = e.get("direction", "both")
        if direction in ("forward", "both"):
            G.add_edge(src_id, tgt_id, edge_type=kind)
            edge_metadata[(src_id, tgt_id)] = meta
        if direction in ("backward", "both"):
            G.add_edge(tgt_id, src_id, edge_type=kind)
            edge_metadata[(tgt_id, src_id)] = meta

    result = {
        "G": G,
        "ground_truth": dict(ground_truth),
        "side_effect_counts": dict(side_effect_counts),
        "edge_metadata": edge_metadata,
        "node_map": node_map,
    }

    print(f"  Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    print(f"  Ground truth: {len(ground_truth)} compounds, "
          f"{sum(len(v) for v in ground_truth.values())} pairs")
    print(f"  Side effects: {len(side_effect_counts)} compounds")

    # Save cache
    print("  Saving enriched cache...")
    with open(enriched_pkl, "wb") as f:
        pickle.dump(result, f, protocol=pickle.HIGHEST_PROTOCOL)

    return result


# ---------------------------------------------------------------------------
# Part 1: Download external mapping data
# ---------------------------------------------------------------------------
def load_drugbank_chembl_map() -> dict[str, list[str]]:
    """DrugBank ID -> list of ChEMBL molecule IDs."""
    path = CACHE_DIR / "drugbank_chembl_map.tsv"
    if not path.exists():
        _download(DRUGBANK_CHEMBL_MAP_URL, path, "DrugBank→ChEMBL mapping")
    mapping = defaultdict(list)
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            mapping[row["drugbank_id"]].append(row["chembl_id"])
    return dict(mapping)


def load_drugbank_pubchem_map() -> dict[str, str]:
    """DrugBank ID -> PubChem CID (string)."""
    path = CACHE_DIR / "drugbank_pubchem_map.tsv"
    if not path.exists():
        _download(DRUGBANK_PUBCHEM_MAP_URL, path, "DrugBank→PubChem mapping")
    mapping = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            mapping[row["drugbank_id"]] = row["pubchem_cid"]
    return mapping


def load_sider_severity() -> dict[str, float]:
    """
    Load SIDER frequency-tier weighted severity per compound.

    Returns PubChem CID -> average severity per side effect.
    Key insight: average severity per SE differentiates between
    "500 rare postmarketing SEs" (safe, e.g. Aspirin) vs
    "10 very frequent SEs" (unsafe, e.g. Thalidomide).
    """
    path = CACHE_DIR / "meddra_freq.tsv.gz"
    if not path.exists():
        _download(SIDER_FREQ_URL, path, "SIDER frequency data")

    # Parse: accumulate severity per compound using frequency tiers
    compound_severity_sum = defaultdict(float)
    compound_se_count = defaultdict(int)

    with gzip.open(path, "rt") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
            stitch_flat = parts[0]  # CID1XXXXXXXX
            freq_desc = parts[4].strip().lower() if len(parts) > 4 else ""
            meddra_type = parts[7]  # PT or LLT
            if meddra_type != "PT":
                continue  # Only Preferred Terms, avoid double-counting

            # Determine severity from frequency tier label
            tier_weight = None
            for tier_name, w in FREQ_TIER_WEIGHTS.items():
                if tier_name in freq_desc:
                    tier_weight = w
                    break

            if tier_weight is None:
                # Fall back to numeric upper_freq if no tier label matched
                try:
                    upper_freq = float(parts[6]) if parts[6] else 0.0
                except ValueError:
                    upper_freq = 0.0
                tier_weight = upper_freq if upper_freq > 0 else 0.3  # moderate default

            # Convert STITCH flat to PubChem CID
            try:
                pubchem_cid = str(int(stitch_flat[3:]) - 100_000_000)
            except (ValueError, IndexError):
                continue

            compound_severity_sum[pubchem_cid] += tier_weight
            compound_se_count[pubchem_cid] += 1

    # Average severity per SE (not total sum — avoids penalizing well-studied drugs)
    result = {}
    for cid, total in compound_severity_sum.items():
        count = compound_se_count[cid]
        result[cid] = total / count  # average tier weight per SE

    print(f"  SIDER severity: {len(result)} compounds "
          f"(avg severity range: {min(result.values()):.3f} - {max(result.values()):.3f})")
    return dict(result)


# ---------------------------------------------------------------------------
# Part 2: ChEMBL bioactivity queries
# ---------------------------------------------------------------------------
def query_chembl_activities(chembl_ids: list[str]) -> dict[str, float]:
    """
    Query ChEMBL REST API for pIC50 data per molecule.
    Returns: chembl_id -> median pchembl_value.
    """
    CHEMBL_CACHE_DIR.mkdir(parents=True, exist_ok=True)

    results = {}
    total = len(chembl_ids)
    t0 = time.time()

    for i, cid in enumerate(chembl_ids):
        cache_file = CHEMBL_CACHE_DIR / f"{cid}.json"

        # Load from cache
        if cache_file.exists():
            with open(cache_file) as f:
                cached = json.load(f)
            if cached.get("pchembl_values"):
                vals = cached["pchembl_values"]
                results[cid] = sorted(vals)[len(vals) // 2]  # median
            continue

        # Query API
        url = CHEMBL_ACTIVITY_URL.format(mol=cid)
        pchembl_values = []

        try:
            page_url = url
            while page_url:
                req = urllib.request.Request(
                    page_url, headers={"Accept": "application/json"}
                )
                with urllib.request.urlopen(req, timeout=30) as resp:
                    data = json.loads(resp.read())

                for act in data.get("activities", []):
                    pval = act.get("pchembl_value")
                    stype = act.get("standard_type", "")
                    srel = act.get("standard_relation", "")
                    if pval and stype in ("IC50", "Ki", "Kd", "EC50"):
                        if srel in ("=", "'='", ""):
                            try:
                                pchembl_values.append(float(pval))
                            except ValueError:
                                pass

                # Pagination
                next_url = data.get("page_meta", {}).get("next")
                if next_url:
                    if next_url.startswith("/"):
                        page_url = "https://www.ebi.ac.uk" + next_url
                    else:
                        page_url = next_url
                else:
                    page_url = None

                time.sleep(0.1)  # Rate limiting

        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
            print(f"    WARN: {cid} query failed: {e}")
            pchembl_values = []

        # Cache result
        with open(cache_file, "w") as f:
            json.dump({"chembl_id": cid, "pchembl_values": pchembl_values}, f)

        if pchembl_values:
            vals = sorted(pchembl_values)
            results[cid] = vals[len(vals) // 2]  # median

        # Progress
        if (i + 1) % 50 == 0 or i == total - 1:
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            remaining = (total - i - 1) / rate if rate > 0 else 0
            print(f"    {i+1}/{total} queried ({len(results)} with data, "
                  f"{elapsed:.0f}s elapsed, ~{remaining:.0f}s remaining)")

    return results


# ---------------------------------------------------------------------------
# Part 3: Build enriched graph weights
# ---------------------------------------------------------------------------
def build_enriched_weights(
    G: nx.DiGraph,
    edge_metadata: dict,
    compound_potency: dict[str, float],  # drugbank_id -> median pIC50
    compound_severity: dict[str, float],  # drugbank_id -> SIDER severity
    side_effect_counts: dict[str, int],   # compound_het_id -> SE count
) -> tuple[nx.DiGraph, nx.DiGraph, nx.DiGraph]:
    """
    Build three weighted graphs (efficacy, safety, evidence) with continuous weights.
    Returns: (G_eff, G_saf, G_evd)
    """
    G_eff = G.copy()
    G_saf = G.copy()
    G_evd = G.copy()

    # Compute normalization constants for efficacy (pIC50)
    all_potencies = list(compound_potency.values())
    if all_potencies:
        sorted_pot = sorted(all_potencies)
        median_potency = sorted_pot[len(sorted_pot) // 2]
        min_potency = sorted_pot[0]
        max_potency = sorted_pot[-1]
        potency_range = max_potency - min_potency
    else:
        median_potency, min_potency, max_potency, potency_range = 6.0, 4.0, 10.0, 6.0

    all_severities = [v for v in compound_severity.values() if v > 0]
    median_severity = sorted(all_severities)[len(all_severities) // 2] if all_severities else 0.5

    all_se_counts = [v for v in side_effect_counts.values() if v > 0]
    median_se = sorted(all_se_counts)[len(all_se_counts) // 2] if all_se_counts else 30.0

    stats = {
        "median_potency": median_potency,
        "potency_range": f"[{min_potency:.2f}, {max_potency:.2f}]",
        "modulation_range": f"±{MODULATION_RANGE*100:.0f}%",
        "median_severity": median_severity,
        "median_se_count": median_se,
        "compounds_with_potency": len(compound_potency),
        "compounds_with_severity": len(compound_severity),
    }

    for u, v in G.edges():
        meta = edge_metadata.get((u, v), {})
        edge_type = meta.get("edge_type", G[u][v].get("edge_type", "unknown"))

        # --- Per-edge metadata signals ---
        pubmed_ids = meta.get("pubmed_ids", [])
        sources = meta.get("sources", [])
        n_pubmed = len(pubmed_ids)
        n_sources = max(len(sources), 1)
        similarity = meta.get("similarity")

        # --- EFFICACY WEIGHT ---
        base_prob = EDGE_TYPE_PROB.get(edge_type, 0.50)

        # For binds edges: use compound pIC50 if available
        if edge_type == "binds":
            compound_id = u if u.startswith("Compound::") else v
            db_id = compound_id.replace("Compound::", "")
            potency = compound_potency.get(db_id)

            if potency is not None:
                # POC v3 fix: pIC50 modulates ±MODULATION_RANGE around base weight
                # Preserves topology-based ranking while adding continuous variation
                modulation = MODULATION_RANGE * (potency - median_potency) / (potency_range + 1e-6)
                p_eff = base_prob * (1.0 + modulation)  # higher pIC50 → higher prob → lower weight
                p_eff = max(0.05, min(0.99, p_eff))
            else:
                # Fallback: use fixed edge-type probability (same as POC1)
                p_eff = base_prob  # 0.80 for binds
        elif edge_type == "resembles" and similarity is not None:
            p_eff = max(0.1, similarity)
        else:
            # Other edge types: fixed type-based probability (matches POC1)
            p_eff = base_prob

        G_eff[u][v]["weight"] = _safe_log(p_eff)

        # --- SAFETY WEIGHT ---
        # Safety is per-compound: find the compound endpoint
        compound_id = None
        for node in (u, v):
            if node.startswith("Compound::"):
                compound_id = node
                break

        if compound_id:
            db_id = compound_id.replace("Compound::", "")
            severity = compound_severity.get(db_id)
            if severity is not None:
                # SIDER severity-weighted
                p_saf = 1.0 / (1.0 + severity / median_severity)
            else:
                # Fallback: Hetionet SE count
                n_se = side_effect_counts.get(compound_id, 0)
                # Use SIDER-calibrated normalization
                p_saf = 1.0 / (1.0 + n_se / median_se)
        else:
            # Non-compound edges: neutral safety
            p_saf = 0.90

        G_saf[u][v]["weight"] = _safe_log(max(p_saf, 0.01))

        # --- EVIDENCE WEIGHT ---
        # Based purely on internal metadata (independent from pIC50)
        if edge_type == "binds":
            # Rich metadata: pubmed + sources + has_action
            has_action = 1.0 if meta.get("actions") else 0.0
            evd_score = (math.log(1 + n_pubmed) + n_sources + has_action) / 8.0
            p_evd = 0.2 + 0.7 * min(evd_score, 1.0)
        elif edge_type == "resembles" and similarity is not None:
            p_evd = max(0.2, similarity)
        elif edge_type == "associates":
            p_evd = 0.3 + 0.3 * min(n_sources / 3.0, 1.0)
        else:
            # Other edges: moderate evidence with small variation
            p_evd = 0.4 + 0.2 * min(n_sources / 3.0, 1.0)

        G_evd[u][v]["weight"] = _safe_log(p_evd)

    return G_eff, G_saf, G_evd, stats


# ---------------------------------------------------------------------------
# Part 4: PharmacotherapyDB loading (same as POC1)
# ---------------------------------------------------------------------------
def load_pharmacotherapydb(
    G: nx.DiGraph, ground_truth: dict[str, list[str]],
) -> tuple[set, set, set]:
    """Returns (positive_pairs, all_compounds, all_diseases).

    Uses ALL 137 Hetionet diseases as evaluation set (not just PharmaDB diseases)
    for fair comparison with POC1.  Compounds are those with ground truth.
    """
    path = CACHE_DIR / "pharmacotherapydb_indications.tsv"
    if not path.exists():
        _download(PHARMACOTHERAPYDB_URL, path, "PharmacotherapyDB")

    positive_pairs = set()

    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row.get("category") != "DM":
                continue
            db_id = row["drugbank_id"]
            doid = row["doid_id"]
            compound = f"Compound::{db_id}"
            disease = f"Disease::{doid}"
            if compound in G and disease in G:
                positive_pairs.add((compound, disease))

    # Use all compounds with ground truth (same as POC1)
    all_compounds = {c for c in ground_truth if c in G}
    # Use ALL Disease nodes in the graph (same as POC1: 137 diseases)
    all_diseases = {n for n in G.nodes() if n.startswith("Disease::")}

    return positive_pairs, all_compounds, all_diseases


# ---------------------------------------------------------------------------
# Part 5: Reverse-graph Dijkstra scoring
# ---------------------------------------------------------------------------
def score_all_pairs(
    G_eff: nx.DiGraph, G_saf: nx.DiGraph, G_evd: nx.DiGraph,
    compounds: set, diseases: set,
) -> dict[tuple[str, str], tuple[float, float, float]]:
    """
    Score all (compound, disease) pairs using reverse-graph Dijkstra.
    Returns: {(compound, disease): (dist_eff, dist_saf, dist_evd)}
    """
    scores = {}
    diseases = sorted(diseases)
    compounds_list = sorted(compounds)

    # Build reverse graphs
    print("    Building reverse graphs...")
    R_eff = G_eff.reverse(copy=True)
    R_saf = G_saf.reverse(copy=True)
    R_evd = G_evd.reverse(copy=True)

    t0 = time.time()
    for i, disease in enumerate(diseases):
        if disease not in R_eff:
            continue

        # Dijkstra from disease on reversed graph = distances TO disease
        try:
            dist_eff = nx.single_source_dijkstra_path_length(R_eff, disease, weight="weight")
        except nx.NetworkXError:
            dist_eff = {}
        try:
            dist_saf = nx.single_source_dijkstra_path_length(R_saf, disease, weight="weight")
        except nx.NetworkXError:
            dist_saf = {}
        try:
            dist_evd = nx.single_source_dijkstra_path_length(R_evd, disease, weight="weight")
        except nx.NetworkXError:
            dist_evd = {}

        for compound in compounds_list:
            d_e = dist_eff.get(compound, float("inf"))
            d_s = dist_saf.get(compound, float("inf"))
            d_v = dist_evd.get(compound, float("inf"))
            scores[(compound, disease)] = (d_e, d_s, d_v)

        if (i + 1) % 20 == 0 or i == 0 or i == len(diseases) - 1:
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed if elapsed > 0 else 1
            remaining = (len(diseases) - i - 1) / rate
            print(f"    Scored {i+1}/{len(diseases)} diseases "
                  f"({elapsed:.0f}s elapsed, ~{remaining:.0f}s remaining)")

    return scores


# ---------------------------------------------------------------------------
# Part 6: AUROC computation
# ---------------------------------------------------------------------------
def compute_auroc(
    scores: dict[tuple[str, str], tuple[float, float, float]],
    positive_pairs: set,
    method: str,
) -> tuple[float, float]:
    """
    Compute AUROC and AUPRC for a given scoring method.
    Returns (auroc, auprc).
    """
    scored_list = []
    for (compound, disease), (d_e, d_s, d_v) in scores.items():
        label = 1 if (compound, disease) in positive_pairs else 0

        # Convert distances to scores (higher = more likely treatment)
        p_e = math.exp(-d_e) if d_e < 50 else 0.0
        p_s = math.exp(-d_s) if d_s < 50 else 0.0
        p_v = math.exp(-d_v) if d_v < 50 else 0.0

        if method == "efficacy_1d":
            score = p_e
        elif method == "eff_safety_2d":
            score = math.exp(-(0.6 * d_e + 0.4 * d_s)) if d_e < 50 and d_s < 50 else 0.0
        elif method == "equal_3d":
            score = math.exp(-(d_e + d_s + d_v) / 3) if all(d < 50 for d in (d_e, d_s, d_v)) else 0.0
        elif method == "weighted_3d":
            score = math.exp(-(0.5 * d_e + 0.3 * d_s + 0.2 * d_v)) if all(d < 50 for d in (d_e, d_s, d_v)) else 0.0
        elif method == "harmonic_3d":
            if p_e > 0 and p_s > 0 and p_v > 0:
                score = 3.0 / (1.0 / p_e + 1.0 / p_s + 1.0 / p_v)
            else:
                score = 0.0
        else:
            score = 0.0

        scored_list.append((score, label))

    # Sort by score descending
    scored_list.sort(key=lambda x: -x[0])

    # AUROC via Wilcoxon-Mann-Whitney
    n_pos = sum(lab for _, lab in scored_list)
    n_neg = len(scored_list) - n_pos
    if n_pos == 0 or n_neg == 0:
        return 0.5, 0.0

    rank_sum = 0.0
    for i, (score, label) in enumerate(scored_list):
        if label == 1:
            rank_sum += (len(scored_list) - i)
    auroc = (rank_sum - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)

    # AUPRC via trapezoidal
    tp = 0
    auprc_sum = 0.0
    for i, (score, label) in enumerate(scored_list):
        if label == 1:
            tp += 1
            precision = tp / (i + 1)
            auprc_sum += precision
    auprc = auprc_sum / n_pos if n_pos > 0 else 0.0

    return auroc, auprc


def bootstrap_ci(
    scores: dict, positive_pairs: set, method: str,
    n_resamples: int = 50, seed: int = 42,
) -> tuple[float, float, float]:
    """Bootstrap 95% CI for AUROC."""
    rng = random.Random(seed)
    pairs = list(scores.keys())
    aurocs = []
    for _ in range(n_resamples):
        sample = rng.choices(pairs, k=len(pairs))
        sample_scores = {p: scores[p] for p in sample}
        sample_pos = positive_pairs & set(sample)
        a, _ = compute_auroc(sample_scores, sample_pos, method)
        aurocs.append(a)
    aurocs.sort()
    lo = aurocs[max(0, int(0.025 * len(aurocs)))]
    hi = aurocs[min(len(aurocs) - 1, int(0.975 * len(aurocs)))]
    return lo, hi


# ---------------------------------------------------------------------------
# Part 7: Pareto analysis
# ---------------------------------------------------------------------------
def compute_pareto_front_2d(points: list[tuple[float, float, int]]) -> list[int]:
    """
    Given [(score1, score2, index), ...] where HIGHER is better,
    return indices of Pareto-optimal points.
    """
    # Sort by score1 descending
    sorted_pts = sorted(points, key=lambda x: -x[0])
    front = []
    max_s2 = -float("inf")
    for s1, s2, idx in sorted_pts:
        if s2 > max_s2:
            front.append(idx)
            max_s2 = s2
    return front


def pareto_analysis(
    scores: dict[tuple[str, str], tuple[float, float, float]],
    positive_pairs: set,
    diseases: set,
    compounds: set,
) -> dict:
    """
    Pareto front analysis per disease.
    Returns metrics and case studies.
    """
    compounds_list = sorted(compounds)
    results = {
        "pareto_recall_at_k": {},  # k -> recall
        "efficacy_recall_at_k": {},
        "pareto_rescue_count": 0,
        "pareto_front_sizes": [],
        "case_studies": [],
    }

    total_pareto_hits = defaultdict(int)
    total_eff_hits = defaultdict(int)
    total_pareto_rescued = 0
    total_positives = 0

    for disease in sorted(diseases):
        # Get scores for all compounds for this disease
        compound_scores = []
        for compound in compounds_list:
            pair = (compound, disease)
            if pair not in scores:
                continue
            d_e, d_s, d_v = scores[pair]
            p_e = math.exp(-d_e) if d_e < 50 else 0.0
            p_s = math.exp(-d_s) if d_s < 50 else 0.0
            is_true = 1 if pair in positive_pairs else 0
            compound_scores.append((compound, p_e, p_s, is_true))

        if not compound_scores:
            continue

        n_pos = sum(c[3] for c in compound_scores)
        if n_pos == 0:
            continue
        total_positives += n_pos

        # Efficacy-only ranking (sorted by p_e descending)
        eff_ranked = sorted(compound_scores, key=lambda x: -x[1])

        # 2D Pareto front (efficacy, safety)
        points = [(c[1], c[2], i) for i, c in enumerate(compound_scores)]
        front_indices = set(compute_pareto_front_2d(points))
        results["pareto_front_sizes"].append(len(front_indices))

        # Pareto ranking: front members first (by efficacy), then non-front by efficacy
        front_compounds = [(compound_scores[i], True) for i in front_indices]
        non_front = [(c, False) for i, c in enumerate(compound_scores) if i not in front_indices]
        front_compounds.sort(key=lambda x: -x[0][1])  # by efficacy
        non_front.sort(key=lambda x: -x[0][1])
        pareto_ranked = front_compounds + non_front

        for k in [10, 20, 50, 100]:
            # Efficacy recall@k
            eff_hits = sum(c[3] for c in eff_ranked[:k])
            total_eff_hits[k] += eff_hits

            # Pareto recall@k
            pareto_hits = sum(c[0][3] for c in pareto_ranked[:k])
            total_pareto_hits[k] += pareto_hits

        # Rescued drugs: in Pareto top-50 but not in efficacy top-50
        eff_top50 = set(c[0] for c in eff_ranked[:50])
        pareto_top50 = set(c[0][0] for c in pareto_ranked[:50])
        for c in pareto_ranked[:50]:
            compound, p_e, p_s, is_true = c[0]
            if is_true and compound not in eff_top50:
                total_pareto_rescued += 1

        # Case study candidates
        if n_pos >= 3 and len(front_indices) >= 5:
            eff_top10_set = set(c[0] for c in eff_ranked[:10])
            for c in pareto_ranked[:20]:
                compound, p_e, p_s, is_true = c[0]
                is_on_front = c[1]
                if is_true and is_on_front and compound not in eff_top10_set:
                    eff_rank = next(
                        (j + 1 for j, ec in enumerate(eff_ranked) if ec[0] == compound),
                        -1
                    )
                    pareto_rank = next(
                        (j + 1 for j, pc in enumerate(pareto_ranked) if pc[0][0] == compound),
                        -1
                    )
                    results["case_studies"].append({
                        "compound": compound,
                        "disease": disease,
                        "eff_rank": eff_rank,
                        "pareto_rank": pareto_rank,
                        "p_efficacy": round(p_e, 6),
                        "p_safety": round(p_s, 6),
                    })

    # Aggregate metrics
    for k in [10, 20, 50, 100]:
        results["pareto_recall_at_k"][k] = total_pareto_hits[k] / total_positives if total_positives else 0
        results["efficacy_recall_at_k"][k] = total_eff_hits[k] / total_positives if total_positives else 0
    results["pareto_rescue_count"] = total_pareto_rescued

    return results


# ---------------------------------------------------------------------------
# Part 8: Dimension correlation analysis
# ---------------------------------------------------------------------------
def analyze_correlations(
    scores: dict[tuple[str, str], tuple[float, float, float]],
) -> dict[str, float]:
    """Compute pairwise Pearson correlations between dimensions."""
    eff_vals, saf_vals, evd_vals = [], [], []
    for (d_e, d_s, d_v) in scores.values():
        if all(d < 50 for d in (d_e, d_s, d_v)):
            eff_vals.append(d_e)
            saf_vals.append(d_s)
            evd_vals.append(d_v)

    def pearson(x, y):
        n = len(x)
        if n < 3:
            return 0.0
        mx = sum(x) / n
        my = sum(y) / n
        cov = sum((xi - mx) * (yi - my) for xi, yi in zip(x, y))
        sx = math.sqrt(sum((xi - mx) ** 2 for xi in x))
        sy = math.sqrt(sum((yi - my) ** 2 for yi in y))
        if sx == 0 or sy == 0:
            return 0.0
        return cov / (sx * sy)

    return {
        "eff_saf": round(pearson(eff_vals, saf_vals), 3),
        "eff_evd": round(pearson(eff_vals, evd_vals), 3),
        "saf_evd": round(pearson(saf_vals, evd_vals), 3),
        "n_valid": len(eff_vals),
    }


# ---------------------------------------------------------------------------
# Part 9: Weight distribution analysis
# ---------------------------------------------------------------------------
def analyze_weight_distributions(
    G_eff: nx.DiGraph, G_saf: nx.DiGraph, G_evd: nx.DiGraph,
    edge_metadata: dict,
) -> dict:
    """Report weight statistics for binds edges to verify continuous variation."""
    stats = {}
    for name, G_dim in [("efficacy", G_eff), ("safety", G_saf), ("evidence", G_evd)]:
        binds_weights = []
        all_weights = []
        for u, v, d in G_dim.edges(data=True):
            w = d.get("weight", 0)
            all_weights.append(w)
            meta = edge_metadata.get((u, v), {})
            if meta.get("edge_type") == "binds":
                binds_weights.append(w)

        if binds_weights:
            binds_weights.sort()
            unique = len(set(round(w, 6) for w in binds_weights))
            stats[name] = {
                "binds_mean": sum(binds_weights) / len(binds_weights),
                "binds_std": (sum((w - sum(binds_weights)/len(binds_weights))**2 for w in binds_weights) / len(binds_weights)) ** 0.5,
                "binds_min": binds_weights[0],
                "binds_max": binds_weights[-1],
                "binds_unique": unique,
                "binds_total": len(binds_weights),
            }
    return stats


# ===========================================================================
# MAIN
# ===========================================================================
def main():
    skip_chembl = "--skip-chembl" in sys.argv

    print("=" * 80)
    print(" ChemPath POC v3: Multi-Objective Pareto Drug Repurposing")
    print(" Fixes: topology-first efficacy, frequency-tier safety")
    print(" Gold standard: PharmacotherapyDB (755 DM indications)")
    print("=" * 80)
    t_start = time.time()

    # ------------------------------------------------------------------
    # Part 0: Load Hetionet with edge metadata
    # ------------------------------------------------------------------
    print("\n[Part 0] Loading Hetionet with edge metadata")
    print("-" * 40)
    het_data = load_hetionet_enriched()
    G = het_data["G"]
    edge_metadata = het_data["edge_metadata"]
    ground_truth = het_data["ground_truth"]
    side_effect_counts = het_data["side_effect_counts"]

    # ------------------------------------------------------------------
    # Part 1: Download external data
    # ------------------------------------------------------------------
    print("\n[Part 1] External data acquisition")
    print("-" * 40)

    # DrugBank -> ChEMBL mapping
    print("  Loading DrugBank→ChEMBL mapping...")
    db_chembl_map = load_drugbank_chembl_map()
    het_compounds = [n for n in G.nodes() if n.startswith("Compound::")]
    het_db_ids = {c.replace("Compound::", "") for c in het_compounds}
    chembl_mapped = sum(1 for db in het_db_ids if db in db_chembl_map)
    print(f"  DrugBank→ChEMBL: {chembl_mapped}/{len(het_db_ids)} compounds ({100*chembl_mapped/len(het_db_ids):.1f}%)")

    # DrugBank -> PubChem mapping (for SIDER)
    print("  Loading DrugBank→PubChem mapping...")
    db_pubchem_map = load_drugbank_pubchem_map()
    pubchem_mapped = sum(1 for db in het_db_ids if db in db_pubchem_map)
    print(f"  DrugBank→PubChem: {pubchem_mapped}/{len(het_db_ids)} compounds ({100*pubchem_mapped/len(het_db_ids):.1f}%)")

    # SIDER 4.1 bias-corrected safety score (Kim's freq × CTCAE-severity)
    print("  Loading SIDER 4.1 safety scores (freq × severity / log-norm)...")
    sider_scores = build_sider_safety_scores()
    fallback_scores = build_fallback_scores(sider_scores)
    all_safety_scores = {**fallback_scores, **sider_scores}  # SIDER takes priority

    # Map Compound::DBXXXXX → DBXXXXX for compatibility with weight builder
    compound_severity = {}
    for hetio_id, score in all_safety_scores.items():
        db_id = hetio_id.replace("Compound::", "")
        if db_id in het_db_ids:
            compound_severity[db_id] = score
    print(f"  SIDER→Hetionet: {len(compound_severity)} compounds with safety data "
          f"(SIDER: {len(sider_scores)}, fallback: {len(fallback_scores)})")

    # ------------------------------------------------------------------
    # Part 2: ChEMBL bioactivity queries
    # ------------------------------------------------------------------
    compound_potency = {}
    if not skip_chembl:
        print("\n[Part 2] ChEMBL bioactivity queries")
        print("-" * 40)

        # Only query ChEMBL IDs for compounds in evaluation set + binds edges
        # (saves ~75% of API calls vs querying all 5K+ ChEMBL IDs)
        eval_db_ids = {c.replace("Compound::", "") for c in ground_truth}
        binds_db_ids = set()
        for u, v, d in G.edges(data=True):
            if d.get("edge_type") == "binds":
                for node in (u, v):
                    if node.startswith("Compound::"):
                        binds_db_ids.add(node.replace("Compound::", ""))
        priority_db_ids = eval_db_ids | binds_db_ids

        all_chembl_ids = set()
        db_to_chembl_list = {}
        for db_id in priority_db_ids:
            chembl_ids = db_chembl_map.get(db_id, [])
            if chembl_ids:
                all_chembl_ids.update(chembl_ids)
                db_to_chembl_list[db_id] = chembl_ids
        print(f"  Unique ChEMBL IDs to query: {len(all_chembl_ids)} "
              f"(from {len(priority_db_ids)} priority compounds)")

        # Query ChEMBL API
        chembl_potency = query_chembl_activities(sorted(all_chembl_ids))
        print(f"  ChEMBL molecules with pIC50 data: {len(chembl_potency)}")

        # Aggregate to DrugBank compound level
        for db_id, chembl_ids in db_to_chembl_list.items():
            pvals = [chembl_potency[cid] for cid in chembl_ids if cid in chembl_potency]
            if pvals:
                compound_potency[db_id] = sorted(pvals)[len(pvals) // 2]

        print(f"  Hetionet compounds with pIC50: {len(compound_potency)}/{len(het_db_ids)}")
    else:
        print("\n[Part 2] ChEMBL queries SKIPPED (--skip-chembl)")
        # Try to load from cache
        cached_count = 0
        eval_db_ids = {c.replace("Compound::", "") for c in ground_truth}
        binds_db_ids = set()
        for u, v, d in G.edges(data=True):
            if d.get("edge_type") == "binds":
                for node in (u, v):
                    if node.startswith("Compound::"):
                        binds_db_ids.add(node.replace("Compound::", ""))
        priority_db_ids = eval_db_ids | binds_db_ids
        for db_id in priority_db_ids:
            chembl_ids = db_chembl_map.get(db_id, [])
            pvals = []
            for cid in chembl_ids:
                cache_file = CHEMBL_CACHE_DIR / f"{cid}.json"
                if cache_file.exists():
                    with open(cache_file) as f:
                        cached = json.load(f)
                    if cached.get("pchembl_values"):
                        vals = cached["pchembl_values"]
                        pvals.append(sorted(vals)[len(vals) // 2])
            if pvals:
                compound_potency[db_id] = sorted(pvals)[len(pvals) // 2]
                cached_count += 1
        print(f"  Loaded {cached_count} compounds from ChEMBL cache")

    # ------------------------------------------------------------------
    # Part 3: Build enriched weights
    # ------------------------------------------------------------------
    print("\n[Part 3] Building enriched graph weights")
    print("-" * 40)
    G_eff, G_saf, G_evd, enrich_stats = build_enriched_weights(
        G, edge_metadata, compound_potency, compound_severity, side_effect_counts
    )
    print(f"  Enrichment stats:")
    for k, v in enrich_stats.items():
        print(f"    {k}: {v}")

    # Weight distribution analysis
    print("\n  Weight distributions (binds edges):")
    wdist = analyze_weight_distributions(G_eff, G_saf, G_evd, edge_metadata)
    for dim, stats in wdist.items():
        print(f"    {dim}: mean={stats['binds_mean']:.3f}, std={stats['binds_std']:.3f}, "
              f"unique={stats['binds_unique']}/{stats['binds_total']}, "
              f"range=[{stats['binds_min']:.3f}, {stats['binds_max']:.3f}]")

    # ------------------------------------------------------------------
    # Part 4: Load PharmacotherapyDB
    # ------------------------------------------------------------------
    print("\n[Part 4] Loading PharmacotherapyDB")
    print("-" * 40)
    positive_pairs, eval_compounds, eval_diseases = load_pharmacotherapydb(G, ground_truth)
    n_pairs = len(eval_compounds) * len(eval_diseases)
    print(f"  Evaluation: {len(eval_compounds)} compounds × {len(eval_diseases)} diseases = {n_pairs:,} pairs")
    print(f"  Positive pairs: {len(positive_pairs)}")
    print(f"  Class ratio: 1:{(n_pairs - len(positive_pairs)) // len(positive_pairs)}")

    # ------------------------------------------------------------------
    # Part 5: Score all pairs
    # ------------------------------------------------------------------
    print("\n[Part 5] Multi-Dimensional Scoring")
    print("-" * 40)
    scores = score_all_pairs(G_eff, G_saf, G_evd, eval_compounds, eval_diseases)
    print(f"    Total scored: {len(scores):,} pairs")

    # ------------------------------------------------------------------
    # Part 6: Dimension correlation analysis
    # ------------------------------------------------------------------
    print("\n[Part 6] Dimension Correlation Analysis")
    print("-" * 40)
    corr = analyze_correlations(scores)
    print(f"  Efficacy ↔ Safety:   r = {corr['eff_saf']:+.3f}")
    print(f"  Efficacy ↔ Evidence: r = {corr['eff_evd']:+.3f}")
    print(f"  Safety ↔ Evidence:   r = {corr['saf_evd']:+.3f}")
    print(f"  (Based on {corr['n_valid']:,} pairs with finite distances in all dims)")

    poc1_corr = {"eff_saf": -0.254, "eff_evd": 0.879, "saf_evd": -0.113}
    print(f"\n  Comparison with POC1 (binary graph):")
    for pair in ["eff_saf", "eff_evd", "saf_evd"]:
        delta = corr[pair] - poc1_corr[pair]
        print(f"    {pair}: {poc1_corr[pair]:+.3f} → {corr[pair]:+.3f} (Δ={delta:+.3f})")

    # ------------------------------------------------------------------
    # Part 7: AUROC Benchmark
    # ------------------------------------------------------------------
    print("\n[Part 7] AUROC Benchmark on PharmacotherapyDB")
    print("-" * 40)

    methods = [
        ("efficacy_1d", "Efficacy only (1D)"),
        ("eff_safety_2d", "Efficacy+Safety (2D)"),
        ("equal_3d", "Equal Geometric (3D)"),
        ("weighted_3d", "Weighted Geometric (3D)"),
        ("harmonic_3d", "Harmonic Mean (3D)"),
    ]

    print(f"\n  {'Method':<30s} | {'AUROC':>17s} | {'AUPRC':>6s} | {'Δ vs 1D':>7s}")
    print(f"  {'─' * 30}─┼─{'─' * 17}─┼─{'─' * 6}─┼─{'─' * 7}")

    baseline_auroc = None
    best_auroc = 0
    best_method = ""

    for method_key, method_name in methods:
        auroc, auprc = compute_auroc(scores, positive_pairs, method_key)
        lo, hi = bootstrap_ci(scores, positive_pairs, method_key)

        if baseline_auroc is None:
            baseline_auroc = auroc
            delta_str = "+0.0000"
        else:
            delta = auroc - baseline_auroc
            delta_str = f"{delta:+.4f}"

        if auroc > best_auroc:
            best_auroc = auroc
            best_method = method_name

        print(f"  {method_name:<30s} | {auroc:.4f}  [{lo:.4f}, {hi:.4f}] | {auprc:.4f} | {delta_str}")

    print(f"\n  Best method: {best_method} (AUROC={best_auroc:.4f})")

    # ------------------------------------------------------------------
    # Part 8: Pareto Analysis
    # ------------------------------------------------------------------
    print("\n[Part 8] Pareto Front Analysis")
    print("-" * 40)
    pareto_results = pareto_analysis(scores, positive_pairs, eval_diseases, eval_compounds)

    print(f"\n  Recall@k comparison (efficacy-only vs Pareto-informed):")
    print(f"  {'k':>5s} | {'Eff-only':>10s} | {'Pareto':>10s} | {'Δ':>8s} | {'Rescued':>8s}")
    print(f"  {'─' * 5}─┼─{'─' * 10}─┼─{'─' * 10}─┼─{'─' * 8}─┼─{'─' * 8}")
    for k in [10, 20, 50, 100]:
        eff_r = pareto_results["efficacy_recall_at_k"].get(k, 0)
        par_r = pareto_results["pareto_recall_at_k"].get(k, 0)
        delta = par_r - eff_r
        print(f"  {k:>5d} | {eff_r:>10.4f} | {par_r:>10.4f} | {delta:>+8.4f} |")

    print(f"\n  Pareto rescue count (top-50): {pareto_results['pareto_rescue_count']} "
          f"true treatments found by Pareto but not by efficacy-only")

    front_sizes = pareto_results["pareto_front_sizes"]
    if front_sizes:
        avg_front = sum(front_sizes) / len(front_sizes)
        print(f"  Avg Pareto front size: {avg_front:.1f} compounds per disease "
              f"(range: {min(front_sizes)}-{max(front_sizes)})")

    # Case studies
    if pareto_results["case_studies"]:
        print(f"\n  Case Studies (drugs rescued by Pareto ranking):")
        print(f"  {'Compound':<25s} | {'Disease':<25s} | {'Eff Rank':>8s} | {'Pareto Rank':>11s} | {'P(eff)':>8s} | {'P(saf)':>8s}")
        print(f"  {'─' * 25}─┼─{'─' * 25}─┼─{'─' * 8}─┼─{'─' * 11}─┼─{'─' * 8}─┼─{'─' * 8}")
        for cs in pareto_results["case_studies"][:10]:
            comp_name = cs["compound"].replace("Compound::", "")
            dis_name = cs["disease"].replace("Disease::", "")
            print(f"  {comp_name:<25s} | {dis_name:<25s} | {cs['eff_rank']:>8d} | {cs['pareto_rank']:>11d} | {cs['p_efficacy']:>8.4f} | {cs['p_safety']:>8.4f}")
    else:
        print("\n  No case studies found (no drugs rescued into top-10 by Pareto)")

    # ------------------------------------------------------------------
    # Verdict
    # ------------------------------------------------------------------
    t_total = time.time() - t_start
    print("\n" + "=" * 80)
    print(" VERDICT")
    print("=" * 80)

    print(f"\n  Baseline (1D efficacy):     AUROC = {baseline_auroc:.4f}")
    print(f"  Best multi-dimensional:     AUROC = {best_auroc:.4f} ({best_method})")
    delta_total = best_auroc - baseline_auroc
    print(f"  Improvement:                Δ = {delta_total:+.4f}")

    # Compare with POC1
    poc1_auroc = 0.7752
    vs_poc1 = best_auroc - poc1_auroc
    print(f"\n  vs POC1 binary graph:       AUROC {poc1_auroc:.4f} → {best_auroc:.4f} (Δ={vs_poc1:+.4f})")

    if delta_total > 0.01:
        print(f"\n  ✓ POSITIVE: Multi-dimensional scoring improves AUROC by {delta_total:+.4f}")
        print(f"    ChEMBL enrichment provides meaningful continuous signal.")
    elif delta_total > -0.01:
        print(f"\n  — NEUTRAL: No meaningful difference ({delta_total:+.4f})")
        if corr["eff_evd"] < 0.5:
            print(f"    Dimensions are now independent (r={corr['eff_evd']:.3f} vs 0.879 in POC1)")
            print(f"    but topology (hop count) still dominates edge weights.")
        else:
            print(f"    Dimensions still correlated (r={corr['eff_evd']:.3f})")
    else:
        print(f"\n  ✗ NEGATIVE: Multi-dimensional scoring hurts ({delta_total:+.4f})")

    pareto_rescue = pareto_results["pareto_rescue_count"]
    if pareto_rescue > 0:
        print(f"\n  Pareto value: {pareto_rescue} true treatments rescued in top-50")
        print(f"    → Multi-objective ranking finds drugs that single-objective misses")
    else:
        print(f"\n  Pareto rescue: 0 drugs rescued")

    print(f"\n  Compounds with ChEMBL pIC50: {len(compound_potency)}")
    print(f"  Compounds with SIDER severity: {len(compound_severity)}")
    print(f"  Total runtime: {t_total:.0f}s")
    print("=" * 80)


if __name__ == "__main__":
    main()
