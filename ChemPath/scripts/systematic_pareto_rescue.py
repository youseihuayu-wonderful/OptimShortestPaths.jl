"""
Systematic Pareto Rescue Analysis

For each of the 137 Hetionet diseases, compute:
1. Single-objective (efficacy-only) top-50 ranking
2. Multi-objective (efficacy + safety) Pareto top-50 ranking
3. Count how many PharmacotherapyDB-validated drugs are "rescued"
   (appear in Pareto top-50 but not in efficacy top-50)

Output: per-disease rescue counts + aggregate statistics for the paper.
"""
from __future__ import annotations

import math
import sys
import time

sys.path.insert(0, __import__("os").path.dirname(__file__))

from chempath_enriched_benchmark import (
    EDGE_TYPE_PROB,
    build_enriched_weights,
    load_hetionet_enriched,
    load_pharmacotherapydb,
    score_all_pairs,
)

try:
    from sider_safety_score import build_sider_safety_scores, build_fallback_scores
except ImportError:
    build_sider_safety_scores = None
    build_fallback_scores = None


def pareto_front_2d(points: list[tuple[str, float, float]]) -> list[tuple[str, float, float]]:
    """Compute 2D Pareto front (minimizing both dimensions).
    points: [(compound_id, efficacy_cost, safety_cost), ...]
    Returns non-dominated points.
    """
    sorted_pts = sorted(points, key=lambda x: x[1])  # sort by efficacy
    front = []
    best_safety = float("inf")
    for p in sorted_pts:
        if p[2] < best_safety:
            front.append(p)
            best_safety = p[2]
    return front


def run_systematic_analysis():
    print("=" * 70)
    print("  SYSTEMATIC PARETO RESCUE ANALYSIS")
    print("  Across all 137 Hetionet diseases")
    print("=" * 70)

    t0 = time.time()
    print("\n[1/4] Loading Hetionet...")
    het_data = load_hetionet_enriched()
    G = het_data["G"]
    edge_metadata = het_data["edge_metadata"]
    side_effect_counts = het_data["side_effect_counts"]

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

    compound_potency = {}  # skip ChEMBL for speed

    ground_truth = het_data["ground_truth"]
    positive_pairs, eval_compounds, eval_diseases = (
        load_pharmacotherapydb(G, ground_truth)
    )
    print(f"  {len(eval_compounds)} eval compounds, {len(eval_diseases)} eval diseases")
    print(f"  {len(positive_pairs)} positive pairs")

    print("\n[3/4] Building enriched weights...")
    G_eff, G_saf, G_evd, stats = build_enriched_weights(
        G, edge_metadata, compound_potency, compound_severity, side_effect_counts
    )

    print("\n[4/4] Scoring all pairs...")
    scores = score_all_pairs(G_eff, G_saf, G_evd, eval_compounds, eval_diseases)
    print(f"  {len(scores)} pairs scored")

    # Analyze per disease
    print("\n" + "=" * 70)
    print("  PER-DISEASE PARETO RESCUE ANALYSIS")
    print("=" * 70)

    TOP_K = 50
    total_rescues = 0
    total_validated_rescues = 0
    diseases_with_rescues = 0
    diseases_with_validated_rescues = 0
    all_rescues = []

    disease_results = []

    for disease in sorted(eval_diseases):
        # Gather all compound scores for this disease
        disease_scores = []
        for (compound, dis), (d_e, d_s, d_v) in scores.items():
            if dis == disease:
                disease_scores.append((compound, d_e, d_s))

        if not disease_scores:
            continue

        # Single-objective top-K (by efficacy distance, ascending)
        eff_ranked = sorted(disease_scores, key=lambda x: x[1])
        eff_top_k = set(c for c, _, _ in eff_ranked[:TOP_K])

        # Multi-objective: Pareto front, then rank by composite
        # First get Pareto front
        front = pareto_front_2d(disease_scores)

        # For all points, compute Pareto-aware ranking:
        # dominated points get penalty, front points get bonus
        front_set = set(c for c, _, _ in front)

        # Composite score: normalize both dimensions, then rank
        if len(disease_scores) > 1:
            max_eff = max(d_e for _, d_e, _ in disease_scores)
            min_eff = min(d_e for _, d_e, _ in disease_scores)
            max_saf = max(d_s for _, _, d_s in disease_scores)
            min_saf = min(d_s for _, _, d_s in disease_scores)
            eff_range = max_eff - min_eff if max_eff > min_eff else 1.0
            saf_range = max_saf - min_saf if max_saf > min_saf else 1.0

            pareto_scored = []
            for c, d_e, d_s in disease_scores:
                norm_e = (d_e - min_eff) / eff_range
                norm_s = (d_s - min_saf) / saf_range
                # 60% efficacy, 40% safety (same as benchmark)
                composite = 0.6 * norm_e + 0.4 * norm_s
                # Pareto front members get small bonus
                if c in front_set:
                    composite -= 0.05
                pareto_scored.append((c, composite))

            pareto_ranked = sorted(pareto_scored, key=lambda x: x[1])
            pareto_top_k = set(c for c, _ in pareto_ranked[:TOP_K])
        else:
            pareto_top_k = eff_top_k

        # Rescued = in Pareto top-K but NOT in efficacy top-K
        rescued = pareto_top_k - eff_top_k

        # Validated rescues = rescued AND in PharmacotherapyDB
        validated_rescued = set()
        for compound in rescued:
            if (compound, disease) in positive_pairs:
                validated_rescued.add(compound)

        n_rescued = len(rescued)
        n_validated = len(validated_rescued)

        if n_rescued > 0:
            diseases_with_rescues += 1
            total_rescues += n_rescued
        if n_validated > 0:
            diseases_with_validated_rescues += 1
            total_validated_rescues += n_validated
            for c in validated_rescued:
                # Find ranks
                eff_rank = next((i+1 for i, (cc, _, _) in enumerate(eff_ranked) if cc == c), "?")
                pareto_rank = next((i+1 for i, (cc, _) in enumerate(pareto_ranked) if cc == c), "?")
                all_rescues.append({
                    "drug": c, "disease": disease,
                    "eff_rank": eff_rank, "pareto_rank": pareto_rank,
                    "delta": eff_rank - pareto_rank if isinstance(eff_rank, int) and isinstance(pareto_rank, int) else "?"
                })

        # Count positive pairs for this disease (for context)
        n_positives = sum(1 for (c, d) in positive_pairs if d == disease)

        disease_results.append({
            "disease": disease,
            "n_compounds": len(disease_scores),
            "n_positives": n_positives,
            "n_rescued": n_rescued,
            "n_validated": n_validated,
        })

    # Print results
    print(f"\n{'Disease':<40s} | {'Cmpds':>5s} | {'Pos':>3s} | {'Rescued':>7s} | {'Valid':>5s}")
    print(f"{'─'*40}─┼─{'─'*5}─┼─{'─'*3}─┼─{'─'*7}─┼─{'─'*5}")
    for r in sorted(disease_results, key=lambda x: -x["n_validated"]):
        if r["n_rescued"] > 0:
            print(f"{r['disease'][:40]:<40s} | {r['n_compounds']:>5d} | {r['n_positives']:>3d} | "
                  f"{r['n_rescued']:>7d} | {r['n_validated']:>5d}")

    # Summary
    print("\n" + "=" * 70)
    print("  AGGREGATE RESULTS")
    print("=" * 70)
    print(f"  Total diseases analyzed:          {len(disease_results)}")
    print(f"  Diseases with any rescue:         {diseases_with_rescues} "
          f"({diseases_with_rescues/len(disease_results)*100:.1f}%)")
    print(f"  Diseases with validated rescue:   {diseases_with_validated_rescues} "
          f"({diseases_with_validated_rescues/len(disease_results)*100:.1f}%)")
    print(f"  Total rescued candidates:         {total_rescues}")
    print(f"  Total validated rescues:          {total_validated_rescues}")
    print(f"  Avg rescues per disease:          {total_rescues/len(disease_results):.1f}")

    if all_rescues:
        print(f"\n  VALIDATED RESCUES (PharmacotherapyDB confirmed):")
        print(f"  {'Drug':<30s} | {'Disease':<30s} | {'Eff→Pareto':>10s} | {'Δ':>5s}")
        print(f"  {'─'*30}─┼─{'─'*30}─┼─{'─'*10}─┼─{'─'*5}")
        for r in sorted(all_rescues, key=lambda x: -(x["delta"] if isinstance(x["delta"], int) else 0)):
            drug_short = r["drug"].replace("Compound::", "")
            disease_short = r["disease"].replace("Disease::", "")[:30]
            print(f"  {drug_short:<30s} | {disease_short:<30s} | "
                  f"{r['eff_rank']:>4}→{r['pareto_rank']:<4} | +{r['delta']}")

    print(f"\n  Total time: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    run_systematic_analysis()
