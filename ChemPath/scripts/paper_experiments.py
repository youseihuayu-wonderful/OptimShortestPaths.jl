"""
Paper Experiments Runner — All analyses for the paper in one script.

Runs on REAL DATA ONLY:
  - Hetionet v1.0 (47,031 nodes, 2.25M edges)
  - PharmacotherapyDB (755 disease-modifying indications)
  - SIDER 4.1 (label-derived frequencies)
  - ChEMBL (pIC50 bioactivity, optional)

NO mock data. NO hard-coded results. Every number is computed live.

Outputs:
  1. AUROC + AUPRC for all methods
  2. Hits@K (K=10, 50, 100) for all methods
  3. Path length analysis (positive vs negative distribution)
  4. 5-fold cross-validation AUROC
  5. p_base sensitivity (5 values)
  6. Systematic Pareto rescue across all 137 diseases
"""
from __future__ import annotations

import json
import math
import os
import pickle
import random
import statistics
import sys
import time
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))

from chempath_enriched_benchmark import (
    EDGE_TYPE_PROB,
    build_enriched_weights,
    compute_auroc,
    compute_pareto_front_2d,
    kfold_cv_auroc,
    load_hetionet_enriched,
    load_pharmacotherapydb,
    analyze_path_lengths,
    score_all_pairs,
    bootstrap_ci,
)

try:
    from sider_safety_score import build_sider_safety_scores, build_fallback_scores
except ImportError:
    build_sider_safety_scores = None
    build_fallback_scores = None

RESULTS_DIR = Path(__file__).resolve().parent.parent / "data" / "paper_results"


def hits_at_k(
    scores: dict[tuple[str, str], tuple[float, float, float]],
    positive_pairs: set,
    method: str,
    k_values: list[int] = [10, 50, 100],
) -> dict[int, float]:
    """Compute Hits@K: fraction of diseases where at least one true positive
    appears in the top-K ranked compounds.

    Uses REAL scores from score_all_pairs, not mock data.
    """
    # Group by disease
    disease_scores: dict[str, list[tuple[str, float, int]]] = {}
    for (compound, disease), (d_e, d_s, d_v) in scores.items():
        p_e = math.exp(-d_e) if d_e < 50 else 0.0
        p_s = math.exp(-d_s) if d_s < 50 else 0.0

        if method == "efficacy_1d":
            score = p_e
        elif method == "eff_safety_2d":
            score = math.exp(-(0.6 * d_e + 0.4 * d_s)) if d_e < 50 and d_s < 50 else 0.0
        else:
            score = p_e

        label = 1 if (compound, disease) in positive_pairs else 0
        if disease not in disease_scores:
            disease_scores[disease] = []
        disease_scores[disease].append((compound, score, label))

    results = {}
    for k in k_values:
        hits = 0
        total_diseases = 0
        for disease, cmpd_scores in disease_scores.items():
            # Only count diseases with at least one positive
            if not any(lab == 1 for _, _, lab in cmpd_scores):
                continue
            total_diseases += 1
            # Sort by score descending
            ranked = sorted(cmpd_scores, key=lambda x: -x[1])[:k]
            if any(lab == 1 for _, _, lab in ranked):
                hits += 1
        results[k] = hits / total_diseases if total_diseases > 0 else 0.0
    return results


def systematic_pareto_rescue(
    scores: dict[tuple[str, str], tuple[float, float, float]],
    positive_pairs: set,
    eval_diseases: set,
    top_k: int = 50,
) -> dict:
    """Systematic Pareto rescue analysis across ALL diseases.

    For each disease:
      - Rank compounds by efficacy (1D)
      - Rank compounds by Pareto-aware composite (2D)
      - Count rescues: in Pareto top-K but not in efficacy top-K
      - Count validated rescues: rescues confirmed in PharmacotherapyDB

    Uses REAL scores, REAL positive_pairs. No mock data.
    """
    total_rescues = 0
    total_validated = 0
    diseases_with_rescues = 0
    diseases_with_validated = 0
    all_validated_rescues = []

    for disease in sorted(eval_diseases):
        disease_scores = []
        for (compound, dis), (d_e, d_s, d_v) in scores.items():
            if dis == disease:
                disease_scores.append((compound, d_e, d_s))

        if len(disease_scores) < 2:
            continue

        # 1D: rank by efficacy distance (ascending = better)
        eff_ranked = sorted(disease_scores, key=lambda x: x[1])
        eff_top_k = set(c for c, _, _ in eff_ranked[:top_k])

        # 2D: Pareto-aware composite
        max_eff = max(d_e for _, d_e, _ in disease_scores)
        min_eff = min(d_e for _, d_e, _ in disease_scores)
        max_saf = max(d_s for _, _, d_s in disease_scores)
        min_saf = min(d_s for _, _, d_s in disease_scores)
        eff_range = max(max_eff - min_eff, 1e-6)
        saf_range = max(max_saf - min_saf, 1e-6)

        # Compute Pareto front for bonus
        # compute_pareto_front_2d expects HIGHER=better, so negate distances
        points = [(-d_e, -d_s, i) for i, (_, d_e, d_s) in enumerate(disease_scores)]
        front_indices = set(compute_pareto_front_2d(points))

        pareto_scored = []
        for i, (c, d_e, d_s) in enumerate(disease_scores):
            norm_e = (d_e - min_eff) / eff_range
            norm_s = (d_s - min_saf) / saf_range
            composite = 0.6 * norm_e + 0.4 * norm_s
            if i in front_indices:
                composite -= 0.05  # Pareto front bonus
            pareto_scored.append((c, composite))

        pareto_ranked = sorted(pareto_scored, key=lambda x: x[1])
        pareto_top_k = set(c for c, _ in pareto_ranked[:top_k])

        rescued = pareto_top_k - eff_top_k
        validated = {c for c in rescued if (c, disease) in positive_pairs}

        if rescued:
            diseases_with_rescues += 1
            total_rescues += len(rescued)
        if validated:
            diseases_with_validated += 1
            total_validated += len(validated)
            for c in validated:
                eff_rank = next((i+1 for i, (cc, _, _) in enumerate(eff_ranked) if cc == c), None)
                par_rank = next((i+1 for i, (cc, _) in enumerate(pareto_ranked) if cc == c), None)
                all_validated_rescues.append({
                    "drug": c, "disease": disease,
                    "eff_rank": eff_rank, "pareto_rank": par_rank,
                    "delta": eff_rank - par_rank if eff_rank and par_rank else 0,
                })

    return {
        "total_diseases": len(eval_diseases),
        "diseases_with_any_rescue": diseases_with_rescues,
        "diseases_with_validated_rescue": diseases_with_validated,
        "total_rescued_candidates": total_rescues,
        "total_validated_rescues": total_validated,
        "validated_rescue_details": sorted(
            all_validated_rescues, key=lambda x: -x["delta"]
        ),
    }


def pbase_sensitivity(
    G, edge_metadata, compound_potency, compound_severity, side_effect_counts,
    positive_pairs, eval_compounds, eval_diseases,
    sweep_values=None,
) -> list[dict]:
    """Sweep p_base(binds) and report AUROC for each value.

    Uses REAL graph data. Only the base probability parameter changes.
    """
    if sweep_values is None:
        sweep_values = [0.50, 0.60, 0.70, 0.80, 0.90]

    results = []
    import copy

    for p_val in sweep_values:
        orig = copy.deepcopy(EDGE_TYPE_PROB)
        EDGE_TYPE_PROB["binds"] = p_val

        G_eff, G_saf, G_evd, _ = build_enriched_weights(
            G, edge_metadata, compound_potency, compound_severity, side_effect_counts
        )
        local_scores = score_all_pairs(G_eff, G_saf, G_evd, eval_compounds, eval_diseases)
        auroc_1d, auprc_1d = compute_auroc(local_scores, positive_pairs, "efficacy_1d")
        auroc_2d, auprc_2d = compute_auroc(local_scores, positive_pairs, "eff_safety_2d")

        for k, v in orig.items():
            EDGE_TYPE_PROB[k] = v

        results.append({
            "p_binds": p_val,
            "auroc_1d": auroc_1d, "auprc_1d": auprc_1d,
            "auroc_2d": auroc_2d, "auprc_2d": auprc_2d,
        })
        print(f"    p_binds={p_val:.2f}  AUROC={auroc_1d:.4f}  AUPRC={auprc_1d:.4f}")

    return results


def main():
    print("=" * 80)
    print("  PAPER EXPERIMENTS — ALL REAL DATA")
    print("  NO mock data. Every number computed live.")
    print("=" * 80)
    t_start = time.time()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    all_results = {}

    # ══════════════════════════════════════════════════════════════
    # LOAD DATA
    # ══════════════════════════════════════════════════════════════
    print("\n[1/8] Loading Hetionet v1.0...")
    het_data = load_hetionet_enriched()
    G = het_data["G"]
    edge_metadata = het_data["edge_metadata"]
    ground_truth = het_data["ground_truth"]
    side_effect_counts = het_data["side_effect_counts"]
    print(f"  Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")

    print("\n[2/8] Loading safety scores (SIDER 4.1 + fallback)...")
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
        print(f"  SIDER: {len(sider_scores)}, Fallback: {len(fallback_scores)}, "
              f"Total: {len(compound_severity)} compounds")
    else:
        print("  WARNING: SIDER not available")

    compound_potency = {}  # Skip ChEMBL for speed (--skip-chembl behavior)

    print("\n[3/8] Loading PharmacotherapyDB gold standard...")
    positive_pairs, eval_compounds, eval_diseases = load_pharmacotherapydb(G, ground_truth)
    print(f"  Eval: {len(eval_compounds)} compounds × {len(eval_diseases)} diseases")
    print(f"  Positive pairs: {len(positive_pairs)} (PharmacotherapyDB DM indications)")

    # ══════════════════════════════════════════════════════════════
    # BUILD WEIGHTS + SCORE (single pass)
    # ══════════════════════════════════════════════════════════════
    print("\n[4/8] Building enriched edge weights...")
    G_eff, G_saf, G_evd, build_stats = build_enriched_weights(
        G, edge_metadata, compound_potency, compound_severity, side_effect_counts
    )

    print("\n[5/8] Scoring all compound-disease pairs...")
    scores = score_all_pairs(G_eff, G_saf, G_evd, eval_compounds, eval_diseases)
    print(f"  Scored {len(scores)} pairs")

    # ══════════════════════════════════════════════════════════════
    # EXPERIMENT 1: AUROC + AUPRC for all methods
    # ══════════════════════════════════════════════════════════════
    print("\n" + "=" * 80)
    print("  EXPERIMENT 1: AUROC + AUPRC (all methods)")
    print("=" * 80)

    methods = ["efficacy_1d", "eff_safety_2d", "weighted_3d", "eff_evd_fused", "eff_evd_safety_fused"]
    method_results = {}
    for m in methods:
        auroc, auprc = compute_auroc(scores, positive_pairs, m)
        lo, hi = bootstrap_ci(scores, positive_pairs, m, n_resamples=500)
        method_results[m] = {"auroc": auroc, "auprc": auprc, "ci_lo": lo, "ci_hi": hi}
        print(f"  {m:<25s}  AUROC={auroc:.4f} [{lo:.4f}, {hi:.4f}]  AUPRC={auprc:.4f}")

    all_results["method_comparison"] = method_results

    # ══════════════════════════════════════════════════════════════
    # EXPERIMENT 2: Hits@K
    # ══════════════════════════════════════════════════════════════
    print("\n" + "=" * 80)
    print("  EXPERIMENT 2: Hits@K (K=10, 50, 100)")
    print("=" * 80)

    hits_results = {}
    for m in ["efficacy_1d", "eff_safety_2d"]:
        hk = hits_at_k(scores, positive_pairs, m, [10, 50, 100])
        hits_results[m] = hk
        print(f"  {m:<25s}  H@10={hk[10]:.3f}  H@50={hk[50]:.3f}  H@100={hk[100]:.3f}")

    all_results["hits_at_k"] = hits_results

    # ══════════════════════════════════════════════════════════════
    # EXPERIMENT 3: Path length analysis
    # ══════════════════════════════════════════════════════════════
    print("\n" + "=" * 80)
    print("  EXPERIMENT 3: Path Length Distribution (positive vs negative)")
    print("=" * 80)

    path_stats = analyze_path_lengths(scores, positive_pairs, G_eff)
    if "error" not in path_stats:
        print(f"  Positive paths: n={path_stats['n_positive']}, "
              f"mean={path_stats['pos_mean']:.3f}, median={path_stats['pos_median']:.3f}")
        print(f"  Negative paths: n={path_stats['n_negative']}, "
              f"mean={path_stats['neg_mean']:.3f}, median={path_stats['neg_median']:.3f}")
        print(f"  Cohen's d = {path_stats['cohens_d']:.3f} ({path_stats['separation']})")
        print(f"  KS statistic = {path_stats['ks_statistic']:.3f}")
    else:
        print(f"  ERROR: {path_stats['error']}")

    all_results["path_length_analysis"] = path_stats

    # ══════════════════════════════════════════════════════════════
    # EXPERIMENT 4: 5-fold cross-validation
    # ══════════════════════════════════════════════════════════════
    print("\n" + "=" * 80)
    print("  EXPERIMENT 4: 5-Fold Cross-Validation (disease split)")
    print("=" * 80)

    cv_results = kfold_cv_auroc(
        G, edge_metadata, compound_potency, compound_severity, side_effect_counts,
        positive_pairs, eval_compounds, eval_diseases, ground_truth, n_folds=5,
    )
    for m, fold_aurocs in cv_results.items():
        mean_a = statistics.mean(fold_aurocs)
        std_a = statistics.stdev(fold_aurocs) if len(fold_aurocs) > 1 else 0
        print(f"  {m:<25s}  mean={mean_a:.4f} ± {std_a:.4f}  folds={[f'{a:.4f}' for a in fold_aurocs]}")

    all_results["cross_validation"] = {
        m: {"mean": statistics.mean(v), "std": statistics.stdev(v) if len(v) > 1 else 0, "folds": v}
        for m, v in cv_results.items()
    }

    # ══════════════════════════════════════════════════════════════
    # EXPERIMENT 5: p_base sensitivity
    # ══════════════════════════════════════════════════════════════
    print("\n" + "=" * 80)
    print("  EXPERIMENT 5: p_base(binds) Sensitivity Analysis")
    print("=" * 80)

    sens_results = pbase_sensitivity(
        G, edge_metadata, compound_potency, compound_severity, side_effect_counts,
        positive_pairs, eval_compounds, eval_diseases,
    )
    aurocs = [r["auroc_1d"] for r in sens_results]
    spread = max(aurocs) - min(aurocs)
    default_auroc = next(r["auroc_1d"] for r in sens_results if r["p_binds"] == 0.80)
    best = max(sens_results, key=lambda r: r["auroc_1d"])

    print(f"\n  Summary:")
    print(f"    AUROC range: [{min(aurocs):.4f}, {max(aurocs):.4f}], spread={spread:.4f}")
    print(f"    Default (0.80): {default_auroc:.4f}")
    print(f"    Best: p_binds={best['p_binds']:.2f}, AUROC={best['auroc_1d']:.4f}")
    if spread < 0.02:
        print(f"    → ROBUST (spread < 0.02)")
    elif spread < 0.05:
        print(f"    → MODERATELY SENSITIVE")
    else:
        print(f"    → SENSITIVE (spread > 0.05)")

    all_results["pbase_sensitivity"] = sens_results

    # ══════════════════════════════════════════════════════════════
    # EXPERIMENT 6: Systematic Pareto rescue
    # ══════════════════════════════════════════════════════════════
    print("\n" + "=" * 80)
    print("  EXPERIMENT 6: Systematic Pareto Rescue (all 137 diseases)")
    print("=" * 80)

    pareto_results = systematic_pareto_rescue(scores, positive_pairs, eval_diseases)
    print(f"  Total diseases: {pareto_results['total_diseases']}")
    print(f"  Diseases with any rescue: {pareto_results['diseases_with_any_rescue']}")
    print(f"  Diseases with validated rescue: {pareto_results['diseases_with_validated_rescue']}")
    print(f"  Total rescued candidates: {pareto_results['total_rescued_candidates']}")
    print(f"  Total validated rescues: {pareto_results['total_validated_rescues']}")

    if pareto_results["validated_rescue_details"]:
        print(f"\n  Top validated rescues (by rank improvement):")
        print(f"  {'Drug':<30s} {'Disease':<30s} {'1D→Pareto':>10s} {'Δ':>5s}")
        print(f"  {'─'*30} {'─'*30} {'─'*10} {'─'*5}")
        for r in pareto_results["validated_rescue_details"][:20]:
            drug = r["drug"].replace("Compound::", "")
            disease = r["disease"].replace("Disease::", "")[:30]
            print(f"  {drug:<30s} {disease:<30s} "
                  f"{r['eff_rank']:>4}→{r['pareto_rank']:<4} +{r['delta']}")

    all_results["pareto_rescue"] = pareto_results

    # ══════════════════════════════════════════════════════════════
    # SAVE ALL RESULTS
    # ══════════════════════════════════════════════════════════════
    output_path = RESULTS_DIR / "paper_experiments.json"
    with open(output_path, "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\n  All results saved to: {output_path}")

    t_total = time.time() - t_start
    print(f"\n{'='*80}")
    print(f"  ALL EXPERIMENTS COMPLETE — {t_total:.0f}s total")
    print(f"{'='*80}")


if __name__ == "__main__":
    main()
