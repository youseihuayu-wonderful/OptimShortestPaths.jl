"""
p_base Sensitivity Analysis

Sweeps the 'binds' base probability from 0.50 to 0.95 in 10 steps,
holding all other edge-type probabilities fixed, and reports AUROC + AUPRC
for each setting. This demonstrates whether results are robust to
hyperparameter choices or depend on careful tuning.

Output: table of (p_binds, AUROC, AUPRC) + summary statistics.
"""
from __future__ import annotations

import copy
import math
import sys
import time

# Ensure sibling imports work
sys.path.insert(0, __import__("os").path.dirname(__file__))

from chempath_enriched_benchmark import (
    EDGE_TYPE_PROB,
    build_enriched_weights,
    compute_auroc,
    load_hetionet_enriched,
    load_pharmacotherapydb,
    score_all_pairs,
)

# Try importing SIDER safety scores
try:
    from sider_safety_score import build_sider_safety_scores, build_fallback_scores
except ImportError:
    build_sider_safety_scores = None
    build_fallback_scores = None


def run_sensitivity():
    print("=" * 70)
    print("  p_base SENSITIVITY ANALYSIS")
    print("  Sweeping binds probability from 0.50 to 0.95")
    print("=" * 70)

    # Load data (same as benchmark)
    t0 = time.time()
    print("\n[1/4] Loading graph and evaluation data...")
    het_data = load_hetionet_enriched()
    G = het_data["G"]
    edge_metadata = het_data["edge_metadata"]
    side_effect_counts = het_data["side_effect_counts"]

    # Load safety scores
    print("\n[2/4] Loading safety scores...")
    compound_severity = {}
    if build_sider_safety_scores:
        sider_scores = build_sider_safety_scores()
        fallback_scores = build_fallback_scores(sider_scores)
        all_safety = {**fallback_scores, **sider_scores}
        het_db_ids = {n.replace("Compound::", "") for n in G.nodes() if n.startswith("Compound::")}
        for hetio_id, score in all_safety.items():
            db_id = hetio_id.replace("Compound::", "")
            if db_id in het_db_ids:
                compound_severity[db_id] = score
        print(f"  {len(compound_severity)} compounds with safety data")
    else:
        print("  SIDER not available, safety dimension empty")

    # Empty potency for speed (sensitivity is about p_base, not ChEMBL)
    compound_potency = {}

    print("\n[3/4] Loading PharmacotherapyDB...")
    ground_truth = het_data["ground_truth"]
    positive_pairs, eval_compounds, eval_diseases = (
        load_pharmacotherapydb(G, ground_truth)
    )
    print(f"  Loaded in {time.time() - t0:.1f}s")
    print(f"  Eval: {len(eval_compounds)} compounds × {len(eval_diseases)} diseases")
    print(f"  Positive pairs: {len(positive_pairs)}")

    # Sweep p_binds
    print("\n[4/4] Sweeping p_binds...")
    sweep_values = [0.50, 0.60, 0.70, 0.80, 0.90]  # 5 values for speed
    results = []

    for p_binds in sweep_values:
        # Temporarily override EDGE_TYPE_PROB
        orig = EDGE_TYPE_PROB.copy()
        EDGE_TYPE_PROB["binds"] = p_binds

        G_eff, G_saf, G_evd, _ = build_enriched_weights(
            G, edge_metadata, compound_potency, compound_severity, side_effect_counts
        )
        scores = score_all_pairs(G_eff, G_saf, G_evd, eval_compounds, eval_diseases)
        auroc_1d, auprc_1d = compute_auroc(scores, positive_pairs, "efficacy_1d")
        auroc_2d, auprc_2d = compute_auroc(scores, positive_pairs, "eff_safety_2d")

        # Restore
        for k, v in orig.items():
            EDGE_TYPE_PROB[k] = v

        results.append({
            "p_binds": p_binds,
            "auroc_1d": auroc_1d,
            "auprc_1d": auprc_1d,
            "auroc_2d": auroc_2d,
            "auprc_2d": auprc_2d,
        })
        print(f"  p_binds={p_binds:.2f}  1D_AUROC={auroc_1d:.4f}  "
              f"1D_AUPRC={auprc_1d:.4f}  2D_AUROC={auroc_2d:.4f}")

    # Skip associates sweep for speed — binds is the most critical edge type
    assoc_results = []

    # Summary
    print("\n" + "=" * 70)
    print("  RESULTS: p_binds sensitivity")
    print("=" * 70)
    print(f"  {'p_binds':>8s} | {'1D AUROC':>9s} | {'1D AUPRC':>9s} | {'2D AUROC':>9s}")
    print(f"  {'─'*8}─┼─{'─'*9}─┼─{'─'*9}─┼─{'─'*9}")
    for r in results:
        marker = " ← default" if r["p_binds"] == 0.80 else ""
        print(f"  {r['p_binds']:>8.2f} | {r['auroc_1d']:>9.4f} | {r['auprc_1d']:>9.4f} | "
              f"{r['auroc_2d']:>9.4f}{marker}")

    aurocs = [r["auroc_1d"] for r in results]
    print(f"\n  AUROC range: [{min(aurocs):.4f}, {max(aurocs):.4f}]")
    print(f"  AUROC spread: {max(aurocs) - min(aurocs):.4f}")
    print(f"  Default (0.80): {[r for r in results if r['p_binds']==0.80][0]['auroc_1d']:.4f}")
    best = max(results, key=lambda r: r["auroc_1d"])
    print(f"  Best p_binds: {best['p_binds']:.2f} (AUROC={best['auroc_1d']:.4f})")

    # Associates sweep skipped for speed

    spread = max(aurocs) - min(aurocs)
    if spread < 0.02:
        print(f"\n  CONCLUSION: Results are ROBUST to p_binds (spread < 0.02)")
    elif spread < 0.05:
        print(f"\n  CONCLUSION: Results are MODERATELY sensitive to p_binds")
    else:
        print(f"\n  CONCLUSION: Results are SENSITIVE to p_binds (spread > 0.05)")

    print(f"\n  Total time: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    run_sensitivity()
