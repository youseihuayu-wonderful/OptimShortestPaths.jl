"""
Path Multiplicity Diagnostic: does path COUNT and METAPATH DIVERSITY
predict therapeutic association better than shortest distance alone?

Hypothesis: Himmelstein's DWPC (degree-weighted path count across metapaths)
outperforms our shortest-distance approach because multiplicity of paths and
diversity of metapath types carry additional predictive signal.

Key insight: Himmelstein got AUROC 0.85 with the SAME binary Hetionet — so
the gap is NOT about continuous weights. It's about what information we
extract from the graph.

Usage:
    python path_multiplicity_diagnostic.py
"""
from __future__ import annotations

import math
import random
import sys
import time
from collections import defaultdict, deque
from pathlib import Path

import networkx as nx

# Import graph loading utilities from the benchmark script
sys.path.insert(0, str(Path(__file__).parent))
from chempath_enriched_benchmark import load_hetionet_enriched, load_pharmacotherapydb

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
SAMPLE_DISEASES = 20
SEED = 42


# ---------------------------------------------------------------------------
# AUROC (Wilcoxon-Mann-Whitney)
# ---------------------------------------------------------------------------
def _auroc(labels_scores: list[tuple[int, float]]) -> float:
    n_pos = sum(lab for lab, _ in labels_scores)
    n_neg = len(labels_scores) - n_pos
    if n_pos == 0 or n_neg == 0:
        return 0.5
    ranked = sorted(labels_scores, key=lambda x: -x[1])
    rank_sum = 0.0
    for i, (label, _score) in enumerate(ranked):
        if label == 1:
            rank_sum += (len(ranked) - i)
    return (rank_sum - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)


# ---------------------------------------------------------------------------
# Fast BFS-based path counting (no path enumeration!)
# ---------------------------------------------------------------------------
def bfs_count_shortest_paths(G: nx.DiGraph, source: str) -> tuple[dict, dict]:
    """BFS from source. Returns (distances, path_counts) dicts.

    path_counts[v] = number of shortest paths from source to v.
    This is O(V+E) — same as a single BFS, no path enumeration.
    """
    dist = {source: 0}
    count = {source: 1}
    queue = deque([source])
    while queue:
        u = queue.popleft()
        d = dist[u]
        for v in G.neighbors(u):
            if v not in dist:
                dist[v] = d + 1
                count[v] = count[u]
                queue.append(v)
            elif dist[v] == d + 1:
                count[v] += count[u]
    return dist, count


def count_metapath_types_bfs(
    G: nx.DiGraph, source: str, targets: set, max_depth: int = 4,
) -> dict[str, set]:
    """BFS from source tracking metapath signatures (edge-type sequences).

    Returns: {target: set of metapath tuples reaching target at shortest dist}

    This is a breadth-first traversal that tracks the SET of edge-type
    sequences, not individual paths. Much faster than path enumeration
    because we merge states with the same (node, depth, metapath_signature).
    But for practical purposes, we track per-node the set of signatures
    that reach it at shortest distance.
    """
    # state: (node) -> set of edge_type tuples that reached it at shortest dist
    # We track the edge types along the path as a tuple
    node_sigs: dict[str, set[tuple]] = {source: {()}}
    node_dist: dict[str, int] = {source: 0}
    queue = deque([source])

    target_sigs: dict[str, set] = {}

    while queue:
        u = queue.popleft()
        d = node_dist[u]
        if d >= max_depth:
            continue

        u_sigs = node_sigs[u]

        for v in G.neighbors(u):
            edge_type = G.edges[u, v].get("edge_type", "?")
            # Extend each signature with this edge type
            new_sigs = set()
            for sig in u_sigs:
                new_sigs.add(sig + (edge_type,))

            if v not in node_dist:
                node_dist[v] = d + 1
                node_sigs[v] = new_sigs
                queue.append(v)
                if v in targets:
                    target_sigs[v] = set(new_sigs)
            elif node_dist[v] == d + 1:
                node_sigs[v] |= new_sigs
                if v in targets:
                    target_sigs.setdefault(v, set()).update(new_sigs)

    return target_sigs


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    print("=" * 72, flush=True)
    print("  PATH MULTIPLICITY DIAGNOSTIC", flush=True)
    print("  Does path count / metapath diversity beat shortest distance?", flush=True)
    print("=" * 72, flush=True)

    # 1. Load Hetionet
    print("\n[1] Loading Hetionet...", flush=True)
    data = load_hetionet_enriched()
    G: nx.DiGraph = data["G"]
    ground_truth = data["ground_truth"]
    print(f"    Nodes: {G.number_of_nodes():,}  Edges: {G.number_of_edges():,}", flush=True)

    # 2. Load PharmacotherapyDB
    print("\n[2] Loading PharmacotherapyDB...", flush=True)
    positive_pairs, all_compounds, all_diseases = load_pharmacotherapydb(G, ground_truth)
    print(f"    Positive pairs: {len(positive_pairs)}", flush=True)
    print(f"    Compounds: {len(all_compounds)}  Diseases: {len(all_diseases)}", flush=True)

    # 3. Sample diseases
    rng = random.Random(SEED)
    diseases_list = sorted(all_diseases)
    sampled_diseases = set(rng.sample(diseases_list, min(SAMPLE_DISEASES, len(diseases_list))))
    print(f"\n[3] Sampled {len(sampled_diseases)} diseases", flush=True)

    eval_positives = {p for p in positive_pairs if p[1] in sampled_diseases}
    print(f"    Positive pairs in sample: {len(eval_positives)}", flush=True)

    # 4. Compute features using FAST BFS counting
    print(f"\n[4] Computing features (BFS path counting — no enumeration)...", flush=True)

    # Use reverse graph: BFS from disease finds distances TO disease
    R = G.reverse(copy=False)

    # Features: {(compound, disease): (dist, path_count, metapath_div)}
    pair_features: dict[tuple, tuple] = {}

    t0 = time.time()
    compounds_sorted = sorted(all_compounds)

    for di, disease in enumerate(sorted(sampled_diseases)):
        # Fast BFS from disease on reverse graph
        dist_map, count_map = bfs_count_shortest_paths(R, disease)

        # Metapath diversity (forward graph: compound -> disease)
        # Only compute for reachable compounds
        reachable = {c for c in compounds_sorted if c in dist_map}

        # Use reverse direction for metapath counting too
        mp_sigs = count_metapath_types_bfs(R, disease, reachable, max_depth=4)

        for compound in compounds_sorted:
            if compound in dist_map:
                d = dist_map[compound]
                c = count_map[compound]
                # Metapath signatures (from disease→compound on reverse graph = compound→disease paths)
                sigs = mp_sigs.get(compound, set())
                pair_features[(compound, disease)] = (d, c, len(sigs))
            else:
                pair_features[(compound, disease)] = (float("inf"), 0, 0)

        elapsed = time.time() - t0
        print(f"    Disease {di+1}/{len(sampled_diseases)}: "
              f"{len(reachable)} reachable, {elapsed:.1f}s", flush=True)

    print(f"    Done in {time.time()-t0:.1f}s", flush=True)

    # 5. Summary statistics
    print("\n[5] Feature comparison: positives vs negatives", flush=True)
    print("-" * 60, flush=True)

    pos_dists, neg_dists = [], []
    pos_counts, neg_counts = [], []
    pos_divs, neg_divs = [], []

    all_pairs = list(pair_features.keys())

    for pair, (d, c, m) in pair_features.items():
        is_pos = pair in eval_positives
        if d < float("inf"):
            if is_pos:
                pos_dists.append(d)
                pos_counts.append(c)
                pos_divs.append(m)
            else:
                neg_dists.append(d)
                neg_counts.append(c)
                neg_divs.append(m)

    def _fmt(vals):
        if not vals:
            return "no data"
        return f"mean={sum(vals)/len(vals):.2f}, median={sorted(vals)[len(vals)//2]}"

    print(f"    {'Feature':<25s} | {'Positives':>30s} | {'Negatives':>30s}", flush=True)
    print(f"    {'─'*25}─┼─{'─'*30}─┼─{'─'*30}", flush=True)
    print(f"    {'Shortest distance':<25s} | {_fmt(pos_dists):>30s} | {_fmt(neg_dists):>30s}", flush=True)
    print(f"    {'Path count':<25s} | {_fmt(pos_counts):>30s} | {_fmt(neg_counts):>30s}", flush=True)
    print(f"    {'Metapath diversity':<25s} | {_fmt(pos_divs):>30s} | {_fmt(neg_divs):>30s}", flush=True)

    # Effect sizes
    def _cohens_d(pos, neg):
        if not pos or not neg:
            return 0.0
        mp = sum(pos) / len(pos)
        mn = sum(neg) / len(neg)
        vp = sum((x-mp)**2 for x in pos) / len(pos) if len(pos) > 1 else 0
        vn = sum((x-mn)**2 for x in neg) / len(neg) if len(neg) > 1 else 0
        pooled = math.sqrt((vp + vn) / 2) if (vp + vn) > 0 else 1
        return (mp - mn) / pooled

    d_dist = _cohens_d([-d for d in pos_dists], [-d for d in neg_dists])  # negate: shorter=better
    d_count = _cohens_d([math.log1p(c) for c in pos_counts], [math.log1p(c) for c in neg_counts])
    d_div = _cohens_d([math.log1p(m) for m in pos_divs], [math.log1p(m) for m in neg_divs])

    print(f"\n    Effect sizes (Cohen's d, positive = predictive):", flush=True)
    print(f"      Shorter distance:     {d_dist:+.3f} {'***' if abs(d_dist) > 0.8 else '**' if abs(d_dist) > 0.5 else '*' if abs(d_dist) > 0.2 else ''}", flush=True)
    print(f"      More paths:           {d_count:+.3f} {'***' if abs(d_count) > 0.8 else '**' if abs(d_count) > 0.5 else '*' if abs(d_count) > 0.2 else ''}", flush=True)
    print(f"      More metapath types:  {d_div:+.3f} {'***' if abs(d_div) > 0.8 else '**' if abs(d_div) > 0.5 else '*' if abs(d_div) > 0.2 else ''}", flush=True)

    # 6. AUROC comparison
    print("\n[6] AUROC comparison", flush=True)
    print("=" * 72, flush=True)

    def _score_dist(feat):
        return -feat[0] if feat[0] < float("inf") else -1e6

    def _score_count(feat):
        return math.log1p(feat[1]) if feat[1] > 0 else -1e6

    def _score_div(feat):
        return math.log1p(feat[2]) if feat[2] > 0 else -1e6

    def _score_dist_count(feat):
        if feat[0] == float("inf") or feat[1] == 0:
            return -1e6
        return -feat[0] + 0.5 * math.log1p(feat[1])

    def _score_dist_div(feat):
        if feat[0] == float("inf") or feat[2] == 0:
            return -1e6
        return -feat[0] + 0.5 * math.log1p(feat[2])

    def _score_dist_count_div(feat):
        if feat[0] == float("inf") or feat[1] == 0:
            return -1e6
        return -feat[0] + 0.3 * math.log1p(feat[1]) + 0.3 * math.log1p(feat[2])

    methods = [
        ("Distance only (baseline)", _score_dist),
        ("Path count only", _score_count),
        ("Metapath diversity only", _score_div),
        ("Distance + path_count", _score_dist_count),
        ("Distance + metapath_div", _score_dist_div),
        ("Distance + count + div", _score_dist_count_div),
    ]

    results = {}
    for name, fn in methods:
        ls = []
        for pair in all_pairs:
            feat = pair_features[pair]
            label = 1 if pair in eval_positives else 0
            ls.append((label, fn(feat)))
        results[name] = _auroc(ls)

    baseline = results["Distance only (baseline)"]
    print(f"\n  {'Method':<35s} | {'AUROC':>7s} | {'Δ':>8s}", flush=True)
    print(f"  {'─'*35}─┼─{'─'*7}─┼─{'─'*8}", flush=True)
    for name, auroc in results.items():
        delta = auroc - baseline
        d_str = f"{delta:+.4f}" if name != "Distance only (baseline)" else "   —"
        marker = " ***" if delta > 0.02 else " **" if delta > 0.01 else " *" if delta > 0.005 else ""
        print(f"  {name:<35s} | {auroc:.4f} | {d_str:>8s}{marker}", flush=True)

    # 7. Interpretation
    print("\n" + "=" * 72, flush=True)
    print("  DEEP ANALYSIS", flush=True)
    print("=" * 72, flush=True)

    best = max(results, key=results.get)
    improvement = results[best] - baseline

    print(f"\n  Baseline (distance only):  AUROC = {baseline:.4f}", flush=True)
    print(f"  Best method:               {best}", flush=True)
    print(f"  Improvement:               {improvement:+.4f}", flush=True)

    # Key question: is the gap STRUCTURAL or about WEIGHTS?
    print(f"\n  KEY QUESTION: Why does Himmelstein get 0.85 with the same graph?", flush=True)
    print(f"  ─────────────────────────────────────────────────────────", flush=True)

    if improvement > 0.02:
        print(f"\n  ANSWER: The gap is STRUCTURAL — our method loses information.", flush=True)
        print(f"  Shortest path = ONE number per pair. Himmelstein uses ~1200 features", flush=True)
        print(f"  (DWPC across all metapath types up to length 4). Path multiplicity", flush=True)
        print(f"  and metapath diversity carry signal that shortest distance discards.", flush=True)
        print(f"\n  Adding more databases would NOT fix this — it would add more edges", flush=True)
        print(f"  but our method would still collapse all paths into one distance.", flush=True)
        print(f"\n  What WOULD fix it:", flush=True)
        print(f"    1. Compute DWPC-like features (count paths per metapath type)", flush=True)
        print(f"    2. Use path count as a second dimension in scoring", flush=True)
        print(f"    3. Feed path features into a classifier (like Himmelstein does)", flush=True)
    elif improvement > 0.005:
        print(f"\n  ANSWER: Mixed — path counting helps a little but the main gap", flush=True)
        print(f"  likely comes from Himmelstein's DEGREE WEIGHTING (DWPC vs DPC).", flush=True)
        print(f"  DWPC downweights paths through hub genes, providing implicit IDF.", flush=True)
        print(f"\n  Adding databases would help IF they provide per-edge confidence", flush=True)
        print(f"  scores that break degeneracy. But the structural limitation remains.", flush=True)
    else:
        print(f"\n  ANSWER: Surprisingly, path counting alone doesn't help much.", flush=True)
        print(f"  The gap may come from:", flush=True)
        print(f"    1. Degree-weighting (DWPC penalizes hub genes)", flush=True)
        print(f"    2. ML integration (logistic regression learns feature interactions)", flush=True)
        print(f"    3. Longer metapaths (length 4 captures more biology)", flush=True)
        print(f"\n  Adding databases: UNCERTAIN — depends on whether new topology", flush=True)
        print(f"  creates discriminative paths that our shortest-path misses.", flush=True)

    # Top metapath signatures
    print(f"\n  Most common metapath types reaching positive pairs:", flush=True)
    mp_counts: dict[tuple, int] = defaultdict(int)
    for pair in eval_positives:
        feat = pair_features.get(pair)
        if feat and feat[0] < float("inf"):
            d = int(feat[0])
            # We don't have the actual signatures stored, but we have the count
            # Let's at least show the distance distribution
            pass

    # Show distance distribution of positives vs negatives
    print(f"\n  Distance distribution:", flush=True)
    for d in range(1, 6):
        np_ = sum(1 for x in pos_dists if x == d)
        nn_ = sum(1 for x in neg_dists if x == d)
        pct_p = 100*np_/len(pos_dists) if pos_dists else 0
        pct_n = 100*nn_/len(neg_dists) if neg_dists else 0
        bar_p = "#" * int(pct_p / 2)
        bar_n = "#" * int(pct_n / 2)
        print(f"    d={d}: pos {np_:>4d} ({pct_p:5.1f}%) {bar_p}", flush=True)
        print(f"         neg {nn_:>4d} ({pct_n:5.1f}%) {bar_n}", flush=True)

    # Does adding more databases help?
    print(f"\n  WOULD MORE DATABASES HELP?", flush=True)
    print(f"  ─────────────────────────", flush=True)
    print(f"  Current degeneracy: 800K+ edges with <4% unique weights", flush=True)

    if improvement > 0.02:
        print(f"  But the real gap is STRUCTURAL (path counting >> distance).", flush=True)
        print(f"  More databases would add edges and potentially new paths,", flush=True)
        print(f"  but without changing our scoring from 'shortest distance'", flush=True)
        print(f"  to 'path count + diversity', the ceiling remains.", flush=True)
        print(f"\n  VERDICT: More databases = marginal help. Method change = real fix.", flush=True)
    else:
        print(f"  And path counting doesn't help much either.", flush=True)
        print(f"  More databases might help by adding LONGER metapaths", flush=True)
        print(f"  (e.g., Drug→Target→Pathway→Gene→Disease via KEGG).", flush=True)
        print(f"  But the primary fix should be DWPC-style scoring.", flush=True)

    print(flush=True)


if __name__ == "__main__":
    main()
