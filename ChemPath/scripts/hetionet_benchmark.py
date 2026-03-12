"""
Large-scale benchmark: ChemPath on Hetionet (47K nodes, 2.25M edges).

Downloads Hetionet v1.0, loads it as a directed graph, holds out known
Compound-treats-Disease edges as ground truth, and benchmarks:

  1. Accuracy: Can multi-hop pathfinding recover held-out indications?
  2. Runtime: Dijkstra SSSP at 47K nodes — baseline for DMY comparison.

Hetionet source: https://github.com/hetio/hetionet
License: CC0 1.0 (public domain)
"""

import bz2
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

sys.path.insert(0, str(Path(__file__).parent.parent))

# ---------------------------------------------------------------------------
# Data download & loading
# ---------------------------------------------------------------------------

HETIONET_JSON_URL = "https://github.com/hetio/hetionet/raw/main/hetnet/json/hetionet-v1.0.json.bz2"

CACHE_DIR = Path(__file__).parent.parent / "data" / "hetionet"


def download_file(url: str, dest: Path) -> Path:
    """Download a file if not already cached."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        print(f"  [cached] {dest.name}")
        return dest
    print(f"  Downloading {url}...")
    urllib.request.urlretrieve(url, dest)
    print(f"  Saved to {dest} ({dest.stat().st_size / 1024:.0f} KB)")
    return dest


def load_hetionet() -> tuple[nx.DiGraph, dict]:
    """
    Load Hetionet v1.0 from bz2-compressed JSON and return (graph, metadata).

    Returns:
        G: Directed graph with all Hetionet edges EXCEPT Compound-treats-Disease
        metadata: dict with ground_truth, node_counts, edge_counts
    """
    bz2_path = download_file(HETIONET_JSON_URL, CACHE_DIR / "hetionet-v1.0.json.bz2")

    print("  Decompressing and parsing JSON...")
    with bz2.open(bz2_path, "rt") as f:
        hetio = json.load(f)

    G = nx.DiGraph()

    # Load nodes
    print("  Loading nodes...")
    node_types = {}
    for node in hetio["nodes"]:
        kind = node["kind"]
        identifier = node["identifier"]
        name = node["name"]
        node_id = f"{kind}::{identifier}"
        G.add_node(node_id,
                   name=name,
                   node_type=kind,
                   raw_id=str(identifier))
        node_types[node_id] = kind

    # Load edges
    print("  Loading edges...")
    ground_truth = defaultdict(list)
    edge_type_counts = defaultdict(int)
    held_out = 0

    # Edge types relevant to drug→disease inference path.
    # Skip ultra-high-volume edges (Gene-participates-{BP,CC,MF} = 800K+)
    # that aren't on the Compound→Gene→Pathway→Disease critical path.
    RELEVANT_EDGES = {
        "binds",         # Compound→Gene (11K, critical)
        "downregulates", # Compound→Gene (131K)
        "upregulates",   # Compound→Gene (124K)
        "interacts",     # Gene↔Gene PPI (147K)
        "regulates",     # Gene→Gene regulation (266K)
        "participates",  # Gene→Pathway (keep only Pathway, not BP/CC/MF)
        "associates",    # Disease↔Gene (13K, critical)
        "localizes",     # Disease→Anatomy (4K)
        "palliates",     # Compound→Disease (390)
        "resembles",     # Compound↔Compound, Disease↔Disease (7K)
        "includes",      # PharmClass→Compound (1K)
        "expresses",     # Anatomy→Gene (526K — large but connects anatomy)
    }
    # Skip: causes (138K, SideEffect noise), presents (3K, Symptom),
    #        covaries (62K, weak signal)

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

        # For "participates", only keep Gene→Pathway (not BP/CC/MF)
        if rel == "participates" and tgt_kind != "Pathway" and src_kind != "Pathway":
            continue

        weight = _edge_type_to_weight(rel)

        if direction == "forward":
            G.add_edge(source_id, target_id, edge_type=rel, weight=weight)
        elif direction == "backward":
            G.add_edge(target_id, source_id, edge_type=rel, weight=weight)
        else:  # "both"
            G.add_edge(source_id, target_id, edge_type=rel, weight=weight)
            G.add_edge(target_id, source_id, edge_type=rel, weight=weight)

    # Node counts by type
    type_counts = defaultdict(int)
    for ntype in node_types.values():
        type_counts[ntype] += 1

    print(f"  Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    print(f"  Held out {held_out} Compound-treats-Disease edges as ground truth")
    print(f"  Ground truth: {len(ground_truth)} compounds with known indications")
    print(f"  Node types: {dict(type_counts)}")

    metadata = {
        "ground_truth": dict(ground_truth),
        "node_counts": dict(type_counts),
        "edge_type_counts": dict(edge_type_counts),
        "held_out_edges": held_out,
    }

    return G, metadata


def _edge_type_to_weight(rel: str) -> float:
    """
    Convert Hetionet edge type to a probability-based weight.

    Uses -log(p) where p reflects edge reliability/strength.
    Direct mechanistic links get higher probability (lower weight).

    Hetionet abbreviated edge types:
      CbG = Compound-binds-Gene, CdG = Compound-downregulates-Gene,
      CuG = Compound-upregulates-Gene, CtD = Compound-treats-Disease,
      CpD = Compound-palliates-Disease, CrC = Compound-resembles-Compound,
      GiG = Gene-interacts-Gene, GcG = Gene-covaries-Gene,
      Gr>G = Gene-regulates-Gene, GpBP/CC/MF = Gene-participates-*,
      GpPW = Gene-participates-Pathway, DaG = Disease-associates-Gene,
      DlA = Disease-localizes-Anatomy, DrD = Disease-resembles-Disease,
      DpS = Disease-presents-Symptom, AeG/AdG/AuG = Anatomy-*-Gene,
      PCiC = PharmacologicClass-includes-Compound,
      CcSE = Compound-causes-SideEffect
    """
    # High confidence mechanistic links
    high_conf = {
        "binds",         # Compound-binds-Gene (direct target interaction)
        "downregulates", # Compound/Anatomy-downregulates-Gene
        "upregulates",   # Compound/Anatomy-upregulates-Gene
        "interacts",     # Gene-interacts-Gene (PPI)
        "participates",  # Gene-participates-{BP,CC,MF,Pathway}
        "associates",    # Disease-associates-Gene
        "localizes",     # Disease-localizes-Anatomy
    }
    # Medium confidence
    med_conf = {
        "palliates",     # Compound-palliates-Disease
        "resembles",     # Compound/Disease-resembles-*
        "regulates",     # Gene-regulates-Gene
        "expresses",     # Anatomy-expresses-Gene
        "includes",      # PharmacologicClass-includes-Compound
    }
    # Lower confidence / indirect
    low_conf = {
        "causes",        # Compound-causes-SideEffect
        "presents",      # Disease-presents-Symptom
        "covaries",      # Gene-covaries-Gene
    }

    if rel in high_conf:
        p = 0.80
    elif rel in med_conf:
        p = 0.50
    elif rel in low_conf:
        p = 0.30
    else:
        p = 0.50  # default

    return -math.log(p)


# ---------------------------------------------------------------------------
# Benchmark
# ---------------------------------------------------------------------------

def run_hetionet_benchmark(
    G: nx.DiGraph,
    ground_truth: dict[str, list[str]],
    max_drugs: int = 100,
    seed: int = 42,
) -> dict:
    """
    Benchmark: For each compound with known indications, run Dijkstra
    to all disease nodes and check if true indications rank highly.
    """
    rng = random.Random(seed)

    # Get all disease nodes
    disease_nodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "Disease"]
    n_diseases = len(disease_nodes)
    disease_set = set(disease_nodes)

    print(f"  Disease nodes: {n_diseases}")

    # Filter to compounds that are in the graph and have ground truth
    eligible = [cid for cid in ground_truth if cid in G]
    if len(eligible) > max_drugs:
        rng.shuffle(eligible)
        eligible = eligible[:max_drugs]

    print(f"  Evaluating {len(eligible)} compounds (of {len(ground_truth)} with ground truth)")

    # Run benchmark
    reciprocal_ranks = []
    recall_at = {1: 0, 3: 0, 5: 0, 10: 0, 20: 0}
    total_true = 0
    dijkstra_times = []

    for i, compound_id in enumerate(eligible):
        true_diseases = set(ground_truth[compound_id]) & disease_set
        if not true_diseases:
            continue

        total_true += len(true_diseases)

        # Run Dijkstra from this compound
        t0 = time.perf_counter()
        try:
            distances = nx.single_source_dijkstra_path_length(
                G, compound_id, weight="weight"
            )
        except nx.NetworkXError:
            continue
        dijkstra_time = time.perf_counter() - t0
        dijkstra_times.append(dijkstra_time)

        # Rank diseases by distance (lower = better)
        disease_dists = []
        for d in disease_nodes:
            dist = distances.get(d, float("inf"))
            disease_dists.append((d, dist))
        disease_dists.sort(key=lambda x: x[1])
        ranked_diseases = [d for d, _ in disease_dists]

        # Find rank of first true disease
        ranks = []
        for true_d in true_diseases:
            if true_d in ranked_diseases:
                ranks.append(ranked_diseases.index(true_d) + 1)

        if ranks:
            best_rank = min(ranks)
            reciprocal_ranks.append(1.0 / best_rank)
            for k in recall_at:
                hits = sum(1 for r in ranks if r <= k)
                recall_at[k] += hits
        else:
            reciprocal_ranks.append(0.0)

        if (i + 1) % 20 == 0:
            print(f"    Progress: {i + 1}/{len(eligible)} compounds evaluated")

    # Compute metrics
    mrr = sum(reciprocal_ranks) / len(reciprocal_ranks) if reciprocal_ranks else 0
    recall = {k: v / total_true if total_true > 0 else 0 for k, v in recall_at.items()}

    # Random baseline MRR
    random_mrr = sum(1.0 / r for r in range(1, n_diseases + 1)) / n_diseases

    # Dijkstra timing
    avg_dijkstra = sum(dijkstra_times) / len(dijkstra_times) if dijkstra_times else 0
    total_dijkstra = sum(dijkstra_times)

    return {
        "mrr": mrr,
        "recall_at_k": recall,
        "n_evaluated": len(reciprocal_ranks),
        "n_diseases": n_diseases,
        "total_true_pairs": total_true,
        "random_mrr": random_mrr,
        "avg_dijkstra_ms": avg_dijkstra * 1000,
        "total_dijkstra_s": total_dijkstra,
        "dijkstra_times": dijkstra_times,
    }


def run_auroc(
    G: nx.DiGraph,
    ground_truth: dict[str, list[str]],
    max_drugs: int = 50,
    seed: int = 42,
) -> dict:
    """Compute AUROC for drug-disease prediction on Hetionet."""
    rng = random.Random(seed)

    disease_nodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "Disease"]
    disease_set = set(disease_nodes)

    eligible = [cid for cid in ground_truth if cid in G]
    if len(eligible) > max_drugs:
        rng.shuffle(eligible)
        eligible = eligible[:max_drugs]

    scores = []
    labels = []

    for compound_id in eligible:
        true_diseases = set(ground_truth[compound_id]) & disease_set

        try:
            distances = nx.single_source_dijkstra_path_length(
                G, compound_id, weight="weight"
            )
        except nx.NetworkXError:
            continue

        for d in disease_nodes:
            dist = distances.get(d, float("inf"))
            # Convert distance to score (lower dist = higher score)
            score = math.exp(-dist) if dist < float("inf") else 0.0
            label = 1 if d in true_diseases else 0
            scores.append(score)
            labels.append(label)

    # Compute AUROC (Wilcoxon-Mann-Whitney)
    n_pos = sum(labels)
    n_neg = len(labels) - n_pos

    if n_pos == 0 or n_neg == 0:
        return {"auroc": 0.0, "n_positive": n_pos, "n_negative": n_neg,
                "n_total": 0, "baseline_aupr": 0.0}

    pairs = sorted(zip(scores, labels), key=lambda x: -x[0])
    tp = 0
    concordant = 0.0
    for score, label in pairs:
        if label == 1:
            tp += 1
        else:
            concordant += tp

    auroc = concordant / (n_pos * n_neg)
    baseline_aupr = n_pos / len(labels)

    return {
        "auroc": auroc,
        "n_positive": n_pos,
        "n_negative": n_neg,
        "n_total": len(labels),
        "baseline_aupr": baseline_aupr,
    }


# ---------------------------------------------------------------------------
# Runtime scaling benchmark
# ---------------------------------------------------------------------------

def benchmark_dijkstra_scaling(G: nx.DiGraph, seed: int = 42) -> list[dict]:
    """
    Measure Dijkstra SSSP runtime at increasing subgraph sizes.
    This establishes the O(m + n log n) baseline that DMY improves upon.
    """
    rng = random.Random(seed)
    all_nodes = list(G.nodes())
    n_total = len(all_nodes)
    m_total = G.number_of_edges()

    # Test at multiple scales
    fractions = [0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 1.00]
    results = []

    for frac in fractions:
        n_sub = max(100, int(n_total * frac))
        if n_sub > n_total:
            n_sub = n_total

        # Sample a connected subgraph via BFS from a random node
        start = rng.choice(all_nodes)
        visited = set()
        queue = [start]
        while queue and len(visited) < n_sub:
            node = queue.pop(0)
            if node in visited:
                continue
            visited.add(node)
            for neighbor in G.successors(node):
                if neighbor not in visited:
                    queue.append(neighbor)
            for neighbor in G.predecessors(node):
                if neighbor not in visited:
                    queue.append(neighbor)

        sub_nodes = list(visited)
        H = G.subgraph(sub_nodes).copy()
        n = H.number_of_nodes()
        m = H.number_of_edges()

        if n < 10:
            continue

        # Pick a source with outgoing edges
        sources = [nd for nd in H.nodes() if H.out_degree(nd) > 0]
        if not sources:
            continue

        # Time multiple Dijkstra runs
        n_runs = min(10, max(1, 1000 // max(1, n)))
        times = []
        for _ in range(n_runs):
            src = rng.choice(sources)
            t0 = time.perf_counter()
            nx.single_source_dijkstra_path_length(H, src, weight="weight")
            times.append(time.perf_counter() - t0)

        avg_ms = (sum(times) / len(times)) * 1000

        # DMY advantage: O(m log^(2/3) n) vs Dijkstra O(m + n log n)
        # DMY saves the log(n) factor on the heap operations, replacing
        # O(n log n) with O(m log^(2/3) n). The advantage appears when
        # the n log n term (heap operations) dominates over the m term
        # (edge scanning), i.e., in graphs where m is not much larger than n.
        log_n = math.log2(max(2, n))
        dijkstra_heap = n * log_n  # heap operations component
        dmy_overhead = m * (log_n ** (2.0 / 3.0))  # DMY total
        dijkstra_total = m + dijkstra_heap  # Dijkstra total

        results.append({
            "fraction": frac,
            "nodes": n,
            "edges": m,
            "avg_dijkstra_ms": avg_ms,
            "n_runs": n_runs,
            "m_over_n": m / n if n > 0 else 0,
        })

        print(f"    n={n:>7,}  m={m:>9,}  m/n={m/n:>5.0f}  Dijkstra={avg_ms:>8.2f}ms")

    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 80)
    print(" LARGE-SCALE BENCHMARK: ChemPath on Hetionet v1.0")
    print(" (47K nodes, 2.25M edges — real biomedical knowledge graph)")
    print("=" * 80)

    # Step 1: Load Hetionet
    print("\n[1] Loading Hetionet v1.0")
    t0 = time.time()
    G, metadata = load_hetionet()
    t_load = time.time() - t0
    print(f"  Load time: {t_load:.1f}s")

    ground_truth = metadata["ground_truth"]

    # Step 2: Ranking benchmark
    print("\n[2] Drug-Disease Ranking Benchmark")
    print("  (Hold-out: Compound-treats-Disease edges removed from graph)")
    t0 = time.time()
    ranking = run_hetionet_benchmark(G, ground_truth, max_drugs=100)
    t_rank = time.time() - t0

    print(f"\n  Results ({ranking['n_evaluated']} compounds, "
          f"{ranking['n_diseases']} diseases):")
    print(f"    MRR:          {ranking['mrr']:.3f}  "
          f"(random = {ranking['random_mrr']:.3f})")
    print(f"    MRR lift:     {ranking['mrr'] / ranking['random_mrr']:.1f}x over random")
    for k in [1, 3, 5, 10, 20]:
        print(f"    Recall@{k:<2d}:    {ranking['recall_at_k'][k]:.1%}")
    print(f"    Avg Dijkstra: {ranking['avg_dijkstra_ms']:.1f}ms per compound")
    print(f"    Total time:   {t_rank:.1f}s")

    # Step 3: AUROC
    print("\n[3] AUROC: Binary Drug-Disease Classification")
    t0 = time.time()
    auroc_result = run_auroc(G, ground_truth, max_drugs=50)
    t_auroc = time.time() - t0

    print(f"    AUROC:        {auroc_result['auroc']:.3f}  (random = 0.500)")
    print(f"    Positive:     {auroc_result['n_positive']}")
    print(f"    Negative:     {auroc_result['n_negative']}")
    print(f"    Total:        {auroc_result['n_total']}")
    print(f"    Compute time: {t_auroc:.1f}s")

    # Step 4: Runtime scaling
    print("\n[4] Runtime Scaling: Dijkstra vs DMY (Theoretical)")
    print("    Dijkstra: O(m + n log n)")
    print("    DMY:      O(m log^(2/3) n)  [STOC 2025]")
    print()
    scaling = benchmark_dijkstra_scaling(G)

    # Step 5: Comparison table
    print(f"\n{'=' * 80}")
    print(" COMPARISON: ChemPath Across Scales")
    print(f"{'=' * 80}")
    print(f"  {'Dataset':<25s} | {'Nodes':>8s} | {'Edges':>8s} | {'MRR':>6s} | {'AUROC':>6s} | {'Dijkstra':>10s}")
    print(f"  {'─' * 75}")
    print(f"  {'Curated (18 drugs)':<25s} | {'51':>8s} | {'107':>8s} | {'0.691':>6s} | {'0.806':>6s} | {'<0.1ms':>10s}")
    print(f"  {'ChEMBL (1553 drugs)':<25s} | {'1,583':>8s} | {'1,739':>8s} | {'0.733':>6s} | {'0.909':>6s} | {'<0.1ms':>10s}")
    print(f"  {'Hetionet (47K nodes)':<25s} | {G.number_of_nodes():>8,} | {G.number_of_edges():>8,} | {ranking['mrr']:>6.3f} | {auroc_result['auroc']:>6.3f} | {ranking['avg_dijkstra_ms']:>8.1f}ms")

    # Runtime analysis
    n = G.number_of_nodes()
    m = G.number_of_edges()
    log_n = math.log2(n)
    avg_degree = m / n

    print(f"\n  Runtime Analysis:")
    print(f"    n = {n:,},  m = {m:,},  avg degree = {avg_degree:.0f}")
    print(f"    Dijkstra per query: {ranking['avg_dijkstra_ms']:.0f}ms (NetworkX Python)")
    print(f"    Total for {ranking['n_evaluated']} queries: {ranking['total_dijkstra_s']:.0f}s")
    print()
    print(f"  DMY Algorithm Context (STOC 2025):")
    print(f"    DMY achieves O(m log^(2/3) n) for directed SSSP.")
    print(f"    Dijkstra with binary heap: O(m + n log n).")
    print(f"    DMY's advantage is strongest on sparse directed graphs")
    print(f"    where m = O(n), giving O(n log^(2/3) n) vs O(n log n).")
    print(f"    At log2(n)={log_n:.1f}: log^(2/3)(n) = {log_n**(2/3):.1f} vs log(n) = {log_n:.1f}")
    print(f"    → DMY saves factor {log_n / log_n**(2/3):.1f}x on the log term.")
    print()
    print(f"  Practical Impact:")
    print(f"    For interactive drug screening (1000s of SSSP queries),")
    print(f"    total query time matters. At Hetionet scale:")
    print(f"      Dijkstra (Python):  {ranking['avg_dijkstra_ms']:.0f}ms × 1000 = {ranking['avg_dijkstra_ms']:.0f}s")
    print(f"      Dijkstra (Julia):   ~{ranking['avg_dijkstra_ms']/50:.0f}ms × 1000 = {ranking['avg_dijkstra_ms']/50:.1f}s  (est. 50x from compiled)")
    print(f"      DMY (Julia):        faster by log^(1/3)(n) = {log_n**(1/3):.1f}x factor")

    # Verdict
    print(f"\n{'=' * 80}")
    print(" VERDICT")
    print(f"{'=' * 80}")

    log_speedup = log_n ** (1.0 / 3.0)
    if ranking["mrr"] > ranking["random_mrr"] * 1.5 and auroc_result["auroc"] > 0.55:
        print(f"  PASS: Multi-hop pathfinding maintains signal on Hetionet.")
        print(f"  MRR = {ranking['mrr']:.3f} ({ranking['mrr'] / ranking['random_mrr']:.1f}x random)")
        print(f"  AUROC = {auroc_result['auroc']:.3f}")
        print(f"  DMY log-term speedup at this scale: {log_speedup:.1f}x")
    elif auroc_result["auroc"] > 0.50:
        print(f"  MARGINAL: Some signal but weak.")
        print(f"  MRR = {ranking['mrr']:.3f}, AUROC = {auroc_result['auroc']:.3f}")
    else:
        print(f"  FAIL: No discriminative signal on Hetionet.")
        print(f"  MRR = {ranking['mrr']:.3f}, AUROC = {auroc_result['auroc']:.3f}")

    print(f"{'=' * 80}")


if __name__ == "__main__":
    main()
