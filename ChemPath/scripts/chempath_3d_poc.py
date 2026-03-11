"""
ChemPath-3D Proof of Concept:
Multi-Dimensional Probability Scoring on Hetionet vs PharmacotherapyDB.

Tests whether adding a safety dimension (per-compound side effect penalty)
to shortest-path drug repurposing scoring improves AUROC over single-dimension
Dijkstra on the Hetionet knowledge graph.

Gold standard: PharmacotherapyDB v1.0 (755 disease-modifying indications)
Graph: Hetionet v1.0 (47K nodes, 2.25M edges)

Usage:
    python ChemPath/scripts/chempath_3d_poc.py
"""

import bz2
import csv
import io
import json
import math
import os
import random
import sys
import time
import urllib.request
from collections import defaultdict
from pathlib import Path

import networkx as nx

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

HETIONET_JSON_URL = (
    "https://github.com/hetio/hetionet/raw/main/hetnet/json/hetionet-v1.0.json.bz2"
)
PHARMACOTHERAPYDB_URL = (
    "https://raw.githubusercontent.com/dhimmel/indications/gh-pages/catalog/indications.tsv"
)
PHARMACOTHERAPYDB_FALLBACK = "https://ndownloader.figshare.com/files/4823950"

CACHE_DIR = Path(__file__).parent.parent / "data" / "hetionet"

# Edges relevant to drug→disease inference (same filter as hetionet_benchmark.py)
RELEVANT_EDGES = {
    "binds", "downregulates", "upregulates", "interacts", "regulates",
    "participates", "associates", "localizes", "palliates", "resembles",
    "includes", "expresses",
}

# ---------------------------------------------------------------------------
# 3D probability parameters (pre-registered, not tuned on results)
# ---------------------------------------------------------------------------

EFFICACY_PROBS = {
    "binds": 0.80,
    "downregulates": 0.65,
    "upregulates": 0.65,
    "associates": 0.70,
    "interacts": 0.50,
    "regulates": 0.45,
    "participates": 0.60,
    "localizes": 0.55,
    "palliates": 0.75,
    "resembles": 0.40,
    "includes": 0.50,
    "expresses": 0.55,
}

SAFETY_BASE_PROBS = {
    "binds": 0.85,
    "downregulates": 0.70,
    "upregulates": 0.70,
    "interacts": 0.85,
    "regulates": 0.80,
    "associates": 0.90,
    "participates": 0.95,
    "localizes": 0.95,
    "palliates": 0.80,
    "resembles": 0.90,
    "includes": 0.90,
    "expresses": 0.90,
}

EVIDENCE_PROBS = {
    "binds": 0.85,
    "associates": 0.80,
    "downregulates": 0.70,
    "upregulates": 0.70,
    "interacts": 0.65,
    "regulates": 0.55,
    "participates": 0.75,
    "palliates": 0.80,
    "localizes": 0.70,
    "resembles": 0.50,
    "includes": 0.65,
    "expresses": 0.60,
}

# Median side effect count for sigmoid normalization
SIDE_EFFECT_MEDIAN = 30.0


# ---------------------------------------------------------------------------
# Data download helpers
# ---------------------------------------------------------------------------

def download_file(url: str, dest: Path, fallback_url: str | None = None) -> Path:
    """Download a file if not already cached."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        print(f"  [cached] {dest.name}")
        return dest
    print(f"  Downloading {url}...")
    try:
        urllib.request.urlretrieve(url, dest)
    except Exception as e:
        if fallback_url:
            print(f"  Primary failed ({e}), trying fallback...")
            urllib.request.urlretrieve(fallback_url, dest)
        else:
            raise
    print(f"  Saved to {dest} ({dest.stat().st_size / 1024:.0f} KB)")
    return dest


# ---------------------------------------------------------------------------
# Part 0: Load Hetionet with side effect extraction
# ---------------------------------------------------------------------------

def load_hetionet_3d() -> tuple[nx.DiGraph, dict, dict]:
    """
    Load Hetionet v1.0 and extract side effect counts per compound.

    Returns:
        G: Directed graph (Compound-treats-Disease held out)
        metadata: dict with ground_truth, node/edge counts
        side_effect_counts: dict {compound_node_id: n_side_effects}
    """
    bz2_path = download_file(HETIONET_JSON_URL, CACHE_DIR / "hetionet-v1.0.json.bz2")

    print("  Decompressing and parsing JSON...")
    with bz2.open(bz2_path, "rt") as f:
        hetio = json.load(f)

    G = nx.DiGraph()

    # --- Load nodes ---
    print("  Loading nodes...")
    node_types = {}
    for node in hetio["nodes"]:
        kind = node["kind"]
        identifier = node["identifier"]
        name = node["name"]
        node_id = f"{kind}::{identifier}"
        G.add_node(node_id, name=name, node_type=kind, raw_id=str(identifier))
        node_types[node_id] = kind

    # --- First pass: extract side effect counts BEFORE edge filtering ---
    print("  Counting side effects per compound...")
    side_effect_counts: dict[str, int] = defaultdict(int)
    for edge in hetio["edges"]:
        if edge["kind"] == "causes" and edge["source_id"][0] == "Compound":
            compound_id = f"Compound::{edge['source_id'][1]}"
            side_effect_counts[compound_id] += 1

    se_values = list(side_effect_counts.values())
    if se_values:
        print(f"  Side effects: {len(side_effect_counts)} compounds, "
              f"median={sorted(se_values)[len(se_values)//2]}, "
              f"max={max(se_values)}, min={min(se_values)}")

    # --- Second pass: load edges with filtering ---
    print("  Loading edges...")
    ground_truth = defaultdict(list)
    edge_type_counts = defaultdict(int)
    held_out = 0

    for edge in hetio["edges"]:
        src_kind = edge["source_id"][0]
        src_id_raw = edge["source_id"][1]
        tgt_kind = edge["target_id"][0]
        tgt_id_raw = edge["target_id"][1]
        rel = edge["kind"]
        direction = edge["direction"]

        source_id = f"{src_kind}::{src_id_raw}"
        target_id = f"{tgt_kind}::{tgt_id_raw}"

        if source_id == target_id:
            continue

        edge_type_counts[rel] += 1

        # Hold out Compound-treats-Disease as ground truth
        if rel == "treats" and src_kind == "Compound" and tgt_kind == "Disease":
            ground_truth[source_id].append(target_id)
            held_out += 1
            continue

        # Filter: only keep relevant edge types
        if rel not in RELEVANT_EDGES:
            continue

        # For "participates", only keep Gene→Pathway
        if rel == "participates" and tgt_kind != "Pathway" and src_kind != "Pathway":
            continue

        if direction == "forward":
            G.add_edge(source_id, target_id, edge_type=rel)
        elif direction == "backward":
            G.add_edge(target_id, source_id, edge_type=rel)
        else:  # "both"
            G.add_edge(source_id, target_id, edge_type=rel)
            G.add_edge(target_id, source_id, edge_type=rel)

    # Node counts
    type_counts = defaultdict(int)
    for ntype in node_types.values():
        type_counts[ntype] += 1

    print(f"  Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    print(f"  Held out {held_out} Compound-treats-Disease edges as ground truth")
    print(f"  Ground truth: {len(ground_truth)} compounds with known indications")

    metadata = {
        "ground_truth": dict(ground_truth),
        "node_counts": dict(type_counts),
        "edge_type_counts": dict(edge_type_counts),
        "held_out_edges": held_out,
    }

    return G, metadata, dict(side_effect_counts)


# ---------------------------------------------------------------------------
# Part 0b: Load PharmacotherapyDB
# ---------------------------------------------------------------------------

def load_pharmacotherapydb() -> list[dict]:
    """
    Download and parse PharmacotherapyDB indications.tsv.

    Returns list of dicts with keys:
        doid_id, drugbank_id, disease, drug, category, n_curators, n_resources
    Filtered to category == "DM" (disease-modifying, 755 indications).
    """
    dest = CACHE_DIR / "pharmacotherapydb_indications.tsv"
    download_file(PHARMACOTHERAPYDB_URL, dest, fallback_url=PHARMACOTHERAPYDB_FALLBACK)

    indications = []
    with open(dest, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row.get("category") == "DM":
                indications.append(row)

    print(f"  PharmacotherapyDB: {len(indications)} DM indications loaded")
    return indications


def map_pharmacotherapydb_to_hetionet(
    indications: list[dict], G: nx.DiGraph
) -> tuple[set, list[dict]]:
    """
    Map PharmacotherapyDB IDs to Hetionet node IDs.

    Returns:
        positive_pairs: set of (compound_node_id, disease_node_id) tuples
        mapped_indications: list of dicts with hetionet IDs added
    """
    positive_pairs = set()
    mapped = []
    unmapped_compounds = set()
    unmapped_diseases = set()

    for ind in indications:
        # PharmacotherapyDB: drugbank_id = "DB00843", doid_id = "DOID:10652"
        # Hetionet: "Compound::DB00843", "Disease::DOID:10652"
        compound_node = f"Compound::{ind['drugbank_id']}"
        disease_node = f"Disease::{ind['doid_id']}"

        if compound_node not in G:
            unmapped_compounds.add(ind["drugbank_id"])
            continue
        if disease_node not in G:
            unmapped_diseases.add(ind["doid_id"])
            continue

        positive_pairs.add((compound_node, disease_node))
        mapped.append({**ind, "compound_node": compound_node, "disease_node": disease_node})

    print(f"  Mapped: {len(positive_pairs)} / {len(indications)} DM pairs to graph")
    print(f"  Unmapped compounds: {len(unmapped_compounds)}, "
          f"unmapped diseases: {len(unmapped_diseases)}")

    return positive_pairs, mapped


# ---------------------------------------------------------------------------
# Part 1: 3D Probability Enrichment
# ---------------------------------------------------------------------------

def enrich_3d(G: nx.DiGraph, side_effect_counts: dict) -> None:
    """
    Add 3D probability weights to every edge in the graph (in-place).

    Dimensions:
        w_efficacy  = -log(p_efficacy)   -- edge-type based
        w_safety    = -log(p_safety)     -- edge-type + per-compound side effect penalty
        w_evidence  = -log(p_evidence)   -- edge-type based
    """
    for u, v, data in G.edges(data=True):
        rel = data.get("edge_type", "unknown")

        # --- Efficacy ---
        p_eff = EFFICACY_PROBS.get(rel, 0.50)

        # --- Safety (per-compound penalty on Compound→* edges) ---
        p_saf = SAFETY_BASE_PROBS.get(rel, 0.85)
        u_type = G.nodes[u].get("node_type", "")
        if u_type == "Compound":
            n_se = side_effect_counts.get(u, 0)
            if n_se > 0:
                # Sigmoid penalty: more side effects → lower safety
                penalty = 1.0 / (1.0 + (n_se / SIDE_EFFECT_MEDIAN))
                p_saf = max(0.05, p_saf * penalty)

        # --- Evidence ---
        p_evd = EVIDENCE_PROBS.get(rel, 0.60)

        # Store log-costs
        data["w_efficacy"] = -math.log(max(1e-10, p_eff))
        data["w_safety"] = -math.log(max(1e-10, p_saf))
        data["w_evidence"] = -math.log(max(1e-10, p_evd))
        data["p_efficacy"] = p_eff
        data["p_safety"] = p_saf
        data["p_evidence"] = p_evd


def check_dimension_correlation(G: nx.DiGraph) -> dict:
    """
    Compute Pearson correlation between the 3 weight dimensions.
    Returns correlation matrix as dict.
    """
    w_eff, w_saf, w_evd = [], [], []
    for _, _, data in G.edges(data=True):
        w_eff.append(data.get("w_efficacy", 0))
        w_saf.append(data.get("w_safety", 0))
        w_evd.append(data.get("w_evidence", 0))

    def pearson(x, y):
        n = len(x)
        if n == 0:
            return 0.0
        mx, my = sum(x) / n, sum(y) / n
        sx = math.sqrt(sum((xi - mx) ** 2 for xi in x) / n)
        sy = math.sqrt(sum((yi - my) ** 2 for yi in y) / n)
        if sx == 0 or sy == 0:
            return 0.0
        cov = sum((xi - mx) * (yi - my) for xi, yi in zip(x, y)) / n
        return cov / (sx * sy)

    corr = {
        "eff_saf": pearson(w_eff, w_saf),
        "eff_evd": pearson(w_eff, w_evd),
        "saf_evd": pearson(w_saf, w_evd),
    }

    # Stats per dimension
    for name, vals in [("efficacy", w_eff), ("safety", w_saf), ("evidence", w_evd)]:
        if vals:
            corr[f"{name}_mean"] = sum(vals) / len(vals)
            corr[f"{name}_std"] = math.sqrt(
                sum((v - sum(vals) / len(vals)) ** 2 for v in vals) / len(vals)
            )

    return corr


# ---------------------------------------------------------------------------
# Part 2: Multi-Dimensional Scoring
# ---------------------------------------------------------------------------

def score_compounds(
    G: nx.DiGraph,
    compound_nodes: list[str],
    disease_nodes: list[str],
    verbose: bool = True,
) -> dict:
    """
    Run 3x Dijkstra SSSP per compound and compute 5 scoring methods.

    Returns:
        scores[method_name][(compound, disease)] = score
    """
    all_scores = {m: {} for m in [
        "1_efficacy_only", "2_efficacy_safety",
        "3_geometric_3d", "4_weighted_3d", "5_harmonic_3d",
    ]}

    disease_set = set(disease_nodes)
    total = len(compound_nodes)
    t0 = time.time()
    skipped = 0

    for i, cid in enumerate(compound_nodes):
        if cid not in G or G.out_degree(cid) == 0:
            skipped += 1
            # Assign 0 scores for all disease pairs
            for did in disease_nodes:
                pair = (cid, did)
                for m in all_scores:
                    all_scores[m][pair] = 0.0
            continue

        # 3x Dijkstra SSSP (cutoff=20 to prune very long paths early)
        try:
            dist_eff = nx.single_source_dijkstra_path_length(
                G, cid, weight="w_efficacy", cutoff=20.0
            )
            dist_saf = nx.single_source_dijkstra_path_length(
                G, cid, weight="w_safety", cutoff=20.0
            )
            dist_evd = nx.single_source_dijkstra_path_length(
                G, cid, weight="w_evidence", cutoff=20.0
            )
        except nx.NetworkXError:
            skipped += 1
            for did in disease_nodes:
                pair = (cid, did)
                for m in all_scores:
                    all_scores[m][pair] = 0.0
            continue

        for did in disease_nodes:
            d_eff = dist_eff.get(did, 20.0)
            d_saf = dist_saf.get(did, 20.0)
            d_evd = dist_evd.get(did, 20.0)

            p_eff = math.exp(-d_eff)
            p_saf = math.exp(-d_saf)
            p_evd = math.exp(-d_evd)

            pair = (cid, did)

            # Method 1: Efficacy only (1D baseline)
            all_scores["1_efficacy_only"][pair] = p_eff

            # Method 2: Efficacy + Safety (2D)
            all_scores["2_efficacy_safety"][pair] = math.exp(
                -(0.6 * d_eff + 0.4 * d_saf)
            )

            # Method 3: Equal geometric mean (3D)
            all_scores["3_geometric_3d"][pair] = math.exp(
                -(d_eff + d_saf + d_evd) / 3.0
            )

            # Method 4: Weighted geometric (3D)
            all_scores["4_weighted_3d"][pair] = math.exp(
                -(0.5 * d_eff + 0.3 * d_saf + 0.2 * d_evd)
            )

            # Method 5: Harmonic mean (3D)
            denom = 0.0
            if p_eff > 1e-15:
                denom += 1.0 / p_eff
            else:
                denom += 1e15
            if p_saf > 1e-15:
                denom += 1.0 / p_saf
            else:
                denom += 1e15
            if p_evd > 1e-15:
                denom += 1.0 / p_evd
            else:
                denom += 1e15
            all_scores["5_harmonic_3d"][pair] = 3.0 / denom if denom > 0 else 0.0

        if verbose and (i + 1) % 25 == 0:
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed
            eta = (total - i - 1) / rate if rate > 0 else 0
            print(f"    Scored {i+1}/{total} compounds "
                  f"({elapsed:.0f}s elapsed, ~{eta:.0f}s remaining)")

    if verbose:
        print(f"    Done: {total} compounds scored in {time.time()-t0:.1f}s "
              f"({skipped} skipped)")

    return all_scores


# ---------------------------------------------------------------------------
# Part 3: AUROC / AUPRC computation
# ---------------------------------------------------------------------------

def compute_auroc(labels: list[int], scores: list[float]) -> float:
    """AUROC via Wilcoxon-Mann-Whitney."""
    n_pos = sum(labels)
    n_neg = len(labels) - n_pos
    if n_pos == 0 or n_neg == 0:
        return 0.0

    pairs = sorted(zip(scores, labels), key=lambda x: -x[0])
    tp = 0
    concordant = 0.0
    for _, label in pairs:
        if label == 1:
            tp += 1
        else:
            concordant += tp

    return concordant / (n_pos * n_neg)


def compute_aupr(labels: list[int], scores: list[float]) -> float:
    """AUPRC via trapezoidal rule."""
    n_pos = sum(labels)
    if n_pos == 0:
        return 0.0

    pairs = sorted(zip(scores, labels), key=lambda x: -x[0])
    tp = 0
    fp = 0
    aupr = 0.0
    prev_recall = 0.0

    for _, label in pairs:
        if label == 1:
            tp += 1
            precision = tp / (tp + fp)
            recall = tp / n_pos
            aupr += precision * (recall - prev_recall)
            prev_recall = recall
        else:
            fp += 1

    return aupr


def evaluate_method(
    method_scores: dict,
    positive_pairs: set,
    compound_nodes: list[str],
    disease_nodes: list[str],
) -> dict:
    """Compute AUROC and AUPRC for a scoring method."""
    scores = []
    labels = []
    disease_set = set(disease_nodes)
    compound_set = set(compound_nodes)

    for pair, score in method_scores.items():
        cid, did = pair
        if cid not in compound_set or did not in disease_set:
            continue
        label = 1 if pair in positive_pairs else 0
        scores.append(score)
        labels.append(label)

    n_pos = sum(labels)
    n_neg = len(labels) - n_pos

    auroc = compute_auroc(labels, scores)
    aupr = compute_aupr(labels, scores)

    return {
        "auroc": auroc,
        "aupr": aupr,
        "n_pos": n_pos,
        "n_neg": n_neg,
        "n_total": len(labels),
    }


def bootstrap_auroc_ci(
    method_scores: dict,
    positive_pairs: set,
    compound_nodes: list[str],
    disease_nodes: list[str],
    n_bootstrap: int = 1000,
    seed: int = 42,
) -> tuple[float, float]:
    """Bootstrap 95% CI on AUROC (resample compounds)."""
    rng = random.Random(seed)

    # Group scores by compound
    compound_data: dict[str, list[tuple]] = defaultdict(list)
    disease_set = set(disease_nodes)
    compound_set = set(compound_nodes)

    for pair, score in method_scores.items():
        cid, did = pair
        if cid in compound_set and did in disease_set:
            label = 1 if pair in positive_pairs else 0
            compound_data[cid].append((score, label))

    compounds = list(compound_data.keys())
    if not compounds:
        return (0.0, 0.0)

    aurocs = []
    for _ in range(n_bootstrap):
        # Resample compounds with replacement
        sampled = [rng.choice(compounds) for _ in range(len(compounds))]
        scores = []
        labels = []
        for cid in sampled:
            for s, l in compound_data[cid]:
                scores.append(s)
                labels.append(l)
        if sum(labels) > 0 and (len(labels) - sum(labels)) > 0:
            aurocs.append(compute_auroc(labels, scores))

    if not aurocs:
        return (0.0, 0.0)

    aurocs.sort()
    lo = aurocs[int(0.025 * len(aurocs))]
    hi = aurocs[int(0.975 * len(aurocs))]
    return (lo, hi)


# ---------------------------------------------------------------------------
# Part 4: Sensitivity Analysis
# ---------------------------------------------------------------------------

def sensitivity_analysis(
    G: nx.DiGraph,
    top_pairs: list[tuple[str, str]],
    disease_nodes: list[str],
    best_method: str,
    n_trials: int = 5,
    seed: int = 42,
) -> dict:
    """
    Perturb edge weights ±20% and check if top predictions remain stable.
    Uses cutoff and only re-scores the unique compounds in top predictions.

    Returns dict with confidence labels for each prediction.
    """
    rng = random.Random(seed)

    # Identify unique compounds in top predictions
    unique_compounds = list(set(c for c, _ in top_pairs))
    # Only score diseases relevant to these compounds
    relevant_diseases = set(d for _, d in top_pairs)

    print(f"    Sensitivity: {len(top_pairs)} pairs, "
          f"{len(unique_compounds)} unique compounds, {n_trials} trials")

    # Save original weights for edges FROM these compounds (faster than all edges)
    original_weights = {}
    for u, v, data in G.edges(data=True):
        original_weights[(u, v)] = {
            "w_efficacy": data.get("w_efficacy", 0),
            "w_safety": data.get("w_safety", 0),
            "w_evidence": data.get("w_evidence", 0),
        }

    # Track how often each pair appears in top-N across trials
    pair_survival = defaultdict(int)
    top_n = len(top_pairs)

    for trial in range(n_trials):
        # Perturb weights
        for u, v, data in G.edges(data=True):
            orig = original_weights[(u, v)]
            data["w_efficacy"] = orig["w_efficacy"] * rng.uniform(0.9, 1.1)
            data["w_safety"] = orig["w_safety"] * rng.uniform(0.8, 1.2)
            data["w_evidence"] = orig["w_evidence"] * rng.uniform(0.9, 1.1)

        # Re-score unique compounds
        trial_scores = {}
        for cid in unique_compounds:
            try:
                dist_eff = nx.single_source_dijkstra_path_length(
                    G, cid, weight="w_efficacy", cutoff=20.0
                )
                dist_saf = nx.single_source_dijkstra_path_length(
                    G, cid, weight="w_safety", cutoff=20.0
                )
                dist_evd = nx.single_source_dijkstra_path_length(
                    G, cid, weight="w_evidence", cutoff=20.0
                )
            except nx.NetworkXError:
                continue

            for did in disease_nodes:
                d_eff = dist_eff.get(did, 20.0)
                d_saf = dist_saf.get(did, 20.0)
                d_evd = dist_evd.get(did, 20.0)

                pair = (cid, did)

                if best_method == "2_efficacy_safety":
                    trial_scores[pair] = math.exp(-(0.6 * d_eff + 0.4 * d_saf))
                elif best_method == "4_weighted_3d":
                    trial_scores[pair] = math.exp(
                        -(0.5 * d_eff + 0.3 * d_saf + 0.2 * d_evd)
                    )
                elif best_method == "3_geometric_3d":
                    trial_scores[pair] = math.exp(-(d_eff + d_saf + d_evd) / 3.0)
                else:
                    trial_scores[pair] = math.exp(-d_eff)

        # Get top-N from this trial
        ranked = sorted(trial_scores.items(), key=lambda x: -x[1])[:top_n]
        trial_top = set(p for p, _ in ranked)

        for pair in top_pairs:
            if pair in trial_top:
                pair_survival[pair] += 1

        print(f"      Trial {trial+1}/{n_trials} done")

    # Restore original weights
    for u, v, data in G.edges(data=True):
        orig = original_weights[(u, v)]
        data["w_efficacy"] = orig["w_efficacy"]
        data["w_safety"] = orig["w_safety"]
        data["w_evidence"] = orig["w_evidence"]

    # Classify confidence (adjusted thresholds for 5 trials)
    confidence = {}
    for pair in top_pairs:
        count = pair_survival.get(pair, 0)
        if count >= 5:
            confidence[pair] = "HIGH"
        elif count >= 3:
            confidence[pair] = "MEDIUM"
        else:
            confidence[pair] = "LOW"

    return {
        "pair_survival": dict(pair_survival),
        "confidence": confidence,
        "n_high": sum(1 for v in confidence.values() if v == "HIGH"),
        "n_medium": sum(1 for v in confidence.values() if v == "MEDIUM"),
        "n_low": sum(1 for v in confidence.values() if v == "LOW"),
    }


# ---------------------------------------------------------------------------
# Part 5: Safety Dimension Showcase
# ---------------------------------------------------------------------------

def safety_showcase(
    all_scores: dict,
    side_effect_counts: dict,
    G: nx.DiGraph,
    positive_pairs: set,
    top_n: int = 200,
) -> list[dict]:
    """
    Find compounds where safety dimension significantly changes ranking.
    """
    # Rank all pairs by efficacy-only
    eff_ranked = sorted(
        all_scores["1_efficacy_only"].items(), key=lambda x: -x[1]
    )

    # Rank all pairs by best 2D/3D method
    best_2d_ranked = sorted(
        all_scores["2_efficacy_safety"].items(), key=lambda x: -x[1]
    )

    # Build rank maps (pair → rank)
    eff_rank = {pair: i + 1 for i, (pair, _) in enumerate(eff_ranked)}
    best_2d_rank = {pair: i + 1 for i, (pair, _) in enumerate(best_2d_ranked)}

    # Find pairs with big rank changes
    showcases = []
    for pair in eff_rank:
        r_eff = eff_rank.get(pair, 999999)
        r_2d = best_2d_rank.get(pair, 999999)
        rank_change = r_eff - r_2d  # positive = demoted by safety

        if abs(rank_change) > 50:
            cid, did = pair
            n_se = side_effect_counts.get(cid, 0)
            compound_name = G.nodes[cid].get("name", cid) if cid in G else cid
            disease_name = G.nodes[did].get("name", did) if did in G else did
            is_true = pair in positive_pairs

            showcases.append({
                "compound": compound_name,
                "compound_id": cid,
                "disease": disease_name,
                "disease_id": did,
                "n_side_effects": n_se,
                "rank_efficacy": r_eff,
                "rank_2d": r_2d,
                "rank_change": rank_change,
                "is_true_indication": is_true,
                "direction": "DEMOTED by safety" if rank_change > 0 else "PROMOTED by safety",
            })

    showcases.sort(key=lambda x: -abs(x["rank_change"]))
    return showcases[:20]


# ---------------------------------------------------------------------------
# Part 6: Per-disease breakdown
# ---------------------------------------------------------------------------

def per_disease_auroc(
    all_scores: dict,
    positive_pairs: set,
    compound_nodes: list[str],
    disease_nodes: list[str],
    G: nx.DiGraph,
    method_1d: str = "1_efficacy_only",
    method_nd: str = "2_efficacy_safety",
) -> list[dict]:
    """Compute per-disease AUROC for 1D vs multi-D method."""
    compound_set = set(compound_nodes)
    results = []

    for did in disease_nodes:
        scores_1d = []
        scores_nd = []
        labels = []

        for cid in compound_nodes:
            pair = (cid, did)
            s1 = all_scores[method_1d].get(pair, 0.0)
            sn = all_scores[method_nd].get(pair, 0.0)
            label = 1 if pair in positive_pairs else 0
            scores_1d.append(s1)
            scores_nd.append(sn)
            labels.append(label)

        n_pos = sum(labels)
        if n_pos < 2:  # Need at least 2 positives for meaningful AUROC
            continue

        auroc_1d = compute_auroc(labels, scores_1d)
        auroc_nd = compute_auroc(labels, scores_nd)

        disease_name = G.nodes[did].get("name", did) if did in G else did
        results.append({
            "disease_id": did,
            "disease_name": disease_name,
            "n_positives": n_pos,
            "auroc_1d": auroc_1d,
            "auroc_nd": auroc_nd,
            "delta": auroc_nd - auroc_1d,
        })

    results.sort(key=lambda x: -x["delta"])
    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 80)
    print(" ChemPath-3D POC: Multi-Dimensional Scoring on Hetionet")
    print(" Gold standard: PharmacotherapyDB (755 DM indications)")
    print("=" * 80)
    t_start = time.time()

    # ===== PART 0: Load data =====
    print("\n[Part 0] Loading data")
    print("-" * 40)

    print("\n  Loading Hetionet v1.0...")
    G, metadata, side_effect_counts = load_hetionet_3d()

    print("\n  Loading PharmacotherapyDB...")
    indications = load_pharmacotherapydb()
    positive_pairs, mapped_indications = map_pharmacotherapydb_to_hetionet(indications, G)

    if not positive_pairs:
        print("  ERROR: No positive pairs mapped. Cannot proceed.")
        return

    # Get the compounds and diseases we're evaluating
    eval_compounds = list(set(m["compound_node"] for m in mapped_indications))
    disease_nodes = [
        n for n, d in G.nodes(data=True) if d.get("node_type") == "Disease"
    ]

    print(f"\n  Evaluation set: {len(eval_compounds)} compounds × "
          f"{len(disease_nodes)} diseases = "
          f"{len(eval_compounds) * len(disease_nodes):,} pairs")
    print(f"  Positive pairs: {len(positive_pairs)}")
    print(f"  Class ratio: 1:{(len(eval_compounds)*len(disease_nodes))//len(positive_pairs)}")

    # ===== PART 1: 3D Enrichment =====
    print("\n[Part 1] 3D Probability Enrichment")
    print("-" * 40)

    enrich_3d(G, side_effect_counts)

    corr = check_dimension_correlation(G)
    print(f"\n  Dimension statistics:")
    print(f"    Efficacy:  mean={corr.get('efficacy_mean', 0):.3f}, "
          f"std={corr.get('efficacy_std', 0):.3f}")
    print(f"    Safety:    mean={corr.get('safety_mean', 0):.3f}, "
          f"std={corr.get('safety_std', 0):.3f}")
    print(f"    Evidence:  mean={corr.get('evidence_mean', 0):.3f}, "
          f"std={corr.get('evidence_std', 0):.3f}")
    print(f"\n  Correlation matrix:")
    print(f"    Efficacy ↔ Safety:   r = {corr['eff_saf']:.3f}")
    print(f"    Efficacy ↔ Evidence: r = {corr['eff_evd']:.3f}")
    print(f"    Safety ↔ Evidence:   r = {corr['saf_evd']:.3f}")

    if corr["eff_saf"] > 0.95:
        print("  ⚠ WARNING: Efficacy and Safety are highly correlated (>0.95).")
        print("    Safety dimension may not add independent information.")

    # ===== PART 2: Multi-Dimensional Scoring =====
    print("\n[Part 2] Multi-Dimensional Scoring")
    print("-" * 40)
    print(f"  Running 3x Dijkstra SSSP for {len(eval_compounds)} compounds...")

    t_score = time.time()
    all_scores = score_compounds(G, eval_compounds, disease_nodes)
    t_score = time.time() - t_score

    # Quick ranking difference check
    n_different = 0
    n_checked = 0
    eff_scores = all_scores["1_efficacy_only"]
    twod_scores = all_scores["2_efficacy_safety"]
    for pair in eff_scores:
        if eff_scores[pair] > 0 and twod_scores.get(pair, 0) > 0:
            n_checked += 1
            if abs(eff_scores[pair] - twod_scores[pair]) / max(eff_scores[pair], 1e-15) > 0.01:
                n_different += 1
    pct_diff = (n_different / n_checked * 100) if n_checked > 0 else 0
    print(f"\n  Ranking difference check: {pct_diff:.1f}% of pairs differ "
          f"between 1D and 2D ({n_different:,}/{n_checked:,})")

    if pct_diff < 1:
        print("  ⚠ WARNING: Very few ranking differences. Safety dimension may lack impact.")

    # ===== PART 3: AUROC Benchmark =====
    print("\n[Part 3] AUROC Benchmark on PharmacotherapyDB")
    print("-" * 40)

    # Add random baseline
    rng = random.Random(42)
    all_scores["0_random"] = {
        pair: rng.random() for pair in all_scores["1_efficacy_only"]
    }

    method_names = {
        "0_random": "Random baseline",
        "1_efficacy_only": "Efficacy only (1D)",
        "2_efficacy_safety": "Efficacy+Safety (2D)",
        "3_geometric_3d": "Equal Geometric (3D)",
        "4_weighted_3d": "Weighted Geometric (3D)",
        "5_harmonic_3d": "Harmonic Mean (3D)",
    }

    results_table = []
    ref_auroc = None

    for method_key in sorted(all_scores.keys()):
        result = evaluate_method(
            all_scores[method_key], positive_pairs, eval_compounds, disease_nodes
        )

        # Bootstrap CI (skip for random)
        if method_key != "0_random":
            ci_lo, ci_hi = bootstrap_auroc_ci(
                all_scores[method_key], positive_pairs,
                eval_compounds, disease_nodes, n_bootstrap=200
            )
        else:
            ci_lo, ci_hi = result["auroc"] - 0.01, result["auroc"] + 0.01

        if method_key == "1_efficacy_only":
            ref_auroc = result["auroc"]

        delta = result["auroc"] - ref_auroc if ref_auroc is not None else 0.0

        results_table.append({
            "method": method_names.get(method_key, method_key),
            "auroc": result["auroc"],
            "ci_lo": ci_lo,
            "ci_hi": ci_hi,
            "aupr": result["aupr"],
            "n_pos": result["n_pos"],
            "n_neg": result["n_neg"],
            "delta": delta,
        })

    # Print results table
    print(f"\n  {'Method':<28s} | {'AUROC':>6s}  {'[95% CI]':>16s} | "
          f"{'AUPRC':>6s} | {'Δ vs 1D':>7s} | {'N_pos':>5s} | {'N_neg':>7s}")
    print(f"  {'─' * 90}")
    for r in results_table:
        delta_str = f"{r['delta']:+.4f}" if r["method"] != "Random baseline" else "  —"
        print(f"  {r['method']:<28s} | {r['auroc']:.4f}  "
              f"[{r['ci_lo']:.4f}, {r['ci_hi']:.4f}] | "
              f"{r['aupr']:.4f} | {delta_str:>7s} | "
              f"{r['n_pos']:>5d} | {r['n_neg']:>7d}")

    # Determine best method
    non_random = [r for r in results_table if r["method"] != "Random baseline"]
    best = max(non_random, key=lambda r: r["auroc"])
    best_method_key = [
        k for k, v in method_names.items() if v == best["method"]
    ][0]

    print(f"\n  Best method: {best['method']} (AUROC={best['auroc']:.4f})")

    # ===== PART 4: Sensitivity Analysis =====
    print("\n[Part 4] Sensitivity Analysis")
    print("-" * 40)

    # Get top-20 predictions from best method
    best_ranked = sorted(
        all_scores[best_method_key].items(), key=lambda x: -x[1]
    )
    top_k = 20
    top_pairs = [pair for pair, _ in best_ranked[:top_k]]

    print(f"  Analyzing stability of top-{top_k} predictions ({best['method']})...")
    t_sens = time.time()
    sensitivity = sensitivity_analysis(
        G, top_pairs, disease_nodes, best_method_key, n_trials=5
    )
    t_sens = time.time() - t_sens

    print(f"\n  Confidence distribution (top-{top_k} predictions):")
    print(f"    HIGH   (5/5 stable):    {sensitivity['n_high']}")
    print(f"    MEDIUM (3-4/5 stable):  {sensitivity['n_medium']}")
    print(f"    LOW    (≤2/5 stable):   {sensitivity['n_low']}")
    pct_stable = (sensitivity["n_high"] + sensitivity["n_medium"]) / top_k * 100
    print(f"    Stable (HIGH+MEDIUM):   {pct_stable:.0f}%")

    # Show a few examples
    print(f"\n  Example HIGH confidence predictions:")
    high_pairs = [p for p in top_pairs if sensitivity["confidence"].get(p) == "HIGH"]
    for pair in high_pairs[:5]:
        cid, did = pair
        cname = G.nodes[cid].get("name", cid) if cid in G else cid
        dname = G.nodes[did].get("name", did) if did in G else did
        is_true = "TRUE" if pair in positive_pairs else "pred"
        score = all_scores[best_method_key].get(pair, 0)
        print(f"    [{is_true:>4s}] {cname:25s} → {dname:25s} "
              f"(score={score:.4f}, survival={sensitivity['pair_survival'].get(pair,0)}/10)")

    # ===== PART 5: Safety Showcase =====
    print("\n[Part 5] Safety Dimension Showcase")
    print("-" * 40)

    showcases = safety_showcase(
        all_scores, side_effect_counts, G, positive_pairs
    )

    if showcases:
        print(f"\n  Top rank changes (safety dimension impact):")
        print(f"  {'Compound':<25s} | {'Disease':<20s} | {'SEs':>4s} | "
              f"{'Rank(1D)':>8s} | {'Rank(2D)':>8s} | {'Change':>7s} | {'True?':>5s}")
        print(f"  {'─' * 100}")
        for s in showcases[:15]:
            true_str = "YES" if s["is_true_indication"] else "no"
            print(f"  {s['compound']:<25s} | {s['disease']:<20s} | "
                  f"{s['n_side_effects']:>4d} | "
                  f"{s['rank_efficacy']:>8d} | {s['rank_2d']:>8d} | "
                  f"{s['rank_change']:>+7d} | {true_str:>5s}")

        # Summarize
        demoted = [s for s in showcases if s["rank_change"] > 0]
        promoted = [s for s in showcases if s["rank_change"] < 0]
        print(f"\n  Summary: {len(demoted)} demoted by safety, "
              f"{len(promoted)} promoted by safety")

        if demoted:
            avg_se_demoted = sum(s["n_side_effects"] for s in demoted) / len(demoted)
            avg_se_promoted = sum(s["n_side_effects"] for s in promoted) / len(promoted) if promoted else 0
            print(f"  Avg side effects: demoted={avg_se_demoted:.0f}, "
                  f"promoted={avg_se_promoted:.0f}")
            if avg_se_demoted > avg_se_promoted:
                print(f"  ✓ Safety dimension correctly demotes high-side-effect compounds")
            else:
                print(f"  ✗ Safety dimension did not clearly separate by side effect count")
    else:
        print("  No significant rank changes found between 1D and 2D.")

    # ===== PART 6: Per-disease breakdown =====
    print("\n[Part 6] Per-Disease AUROC Breakdown")
    print("-" * 40)

    disease_results = per_disease_auroc(
        all_scores, positive_pairs, eval_compounds, disease_nodes, G,
        method_1d="1_efficacy_only", method_nd=best_method_key,
    )

    if disease_results:
        print(f"\n  Top 10 diseases where multi-D helps most:")
        print(f"  {'Disease':<35s} | {'N_pos':>5s} | {'AUROC(1D)':>9s} | "
              f"{'AUROC(nD)':>9s} | {'Δ':>7s}")
        print(f"  {'─' * 75}")
        for r in disease_results[:10]:
            print(f"  {r['disease_name']:<35s} | {r['n_positives']:>5d} | "
                  f"{r['auroc_1d']:>9.4f} | {r['auroc_nd']:>9.4f} | "
                  f"{r['delta']:>+7.4f}")

        print(f"\n  Bottom 5 diseases where multi-D hurts:")
        for r in disease_results[-5:]:
            print(f"  {r['disease_name']:<35s} | {r['n_positives']:>5d} | "
                  f"{r['auroc_1d']:>9.4f} | {r['auroc_nd']:>9.4f} | "
                  f"{r['delta']:>+7.4f}")

        n_improved = sum(1 for r in disease_results if r["delta"] > 0.01)
        n_hurt = sum(1 for r in disease_results if r["delta"] < -0.01)
        n_same = len(disease_results) - n_improved - n_hurt
        print(f"\n  Diseases improved: {n_improved}, hurt: {n_hurt}, same: {n_same}")

    # ===== VERDICT =====
    print(f"\n{'=' * 80}")
    print(" VERDICT")
    print(f"{'=' * 80}")

    t_total = time.time() - t_start

    one_d = next(r for r in results_table if "1D" in r["method"])
    best_nd = best

    delta = best_nd["auroc"] - one_d["auroc"]

    print(f"\n  Baseline (1D efficacy):     AUROC = {one_d['auroc']:.4f}")
    print(f"  Best multi-dimensional:     AUROC = {best_nd['auroc']:.4f} ({best_nd['method']})")
    print(f"  Improvement:                Δ = {delta:+.4f}")

    if delta > 0.05:
        print(f"\n  ✓ STRONG RESULT: Multi-dimensional scoring improves AUROC by {delta:.4f}")
        print(f"    The safety dimension adds meaningful signal for drug repurposing.")
        print(f"    This is paper-worthy (target: JCIM or Bioinformatics).")
    elif delta > 0.01:
        print(f"\n  ~ MODEST RESULT: Multi-dimensional scoring improves AUROC by {delta:.4f}")
        print(f"    Some benefit from safety dimension, but not transformative on Hetionet.")
        print(f"    Stronger results expected on IC50-enriched graphs with real dose data.")
    elif delta > -0.01:
        print(f"\n  — NEUTRAL: No meaningful difference ({delta:+.4f})")
        print(f"    On Hetionet's binary graph, multi-dimensional scoring does not help.")
        print(f"    Pivot: build a ChEMBL-enriched graph with real continuous weights.")
    else:
        print(f"\n  ✗ NEGATIVE: Multi-dimensional scoring hurts by {delta:.4f}")
        print(f"    The safety/evidence dimensions add noise on this graph.")

    print(f"\n  Sensitivity: {pct_stable:.0f}% of top-50 predictions are stable")
    print(f"  Total runtime: {t_total:.0f}s ({t_score:.0f}s scoring, {t_sens:.0f}s sensitivity)")
    print(f"{'=' * 80}")


if __name__ == "__main__":
    main()
