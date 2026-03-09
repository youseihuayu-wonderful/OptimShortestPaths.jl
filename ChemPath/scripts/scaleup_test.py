"""
Scale-up validation: Does multi-hop pathfinding maintain accuracy
when the compound search space grows from 18 to ~200-500?

Fetches real ChEMBL IC50 data for all 10 anticancer targets,
builds the 4-layer graph, and runs the same benchmark.

Key question: If AUROC stays >0.7 with hundreds of distractors,
the framework is viable at larger scale.
"""

import sys
import time
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from chempath.data.curated_data import load_curated_data
from chempath.data.chembl_client import (
    fetch_real_data, enrich_compounds_with_phase,
    save_data, load_saved_data, ANTICANCER_TARGETS,
)
from chempath.graph.network import build_multihop_graph, get_multihop_summary
from chempath.graph.benchmark import (
    run_multihop_benchmark, run_random_baseline, run_single_hop_baseline,
    compute_auroc_aupr, format_benchmark_report, GROUND_TRUTH,
)


DATA_PATH = Path(__file__).parent.parent / "data" / "chembl_scaleup.json"


def fetch_or_load_data(max_per_target: int = 200, use_cache: bool = True) -> dict:
    """Fetch real ChEMBL data, or load from disk if already fetched."""
    if DATA_PATH.exists():
        print(f"  Loading cached data from {DATA_PATH}...")
        return load_saved_data(DATA_PATH)

    print(f"  Fetching IC50 data from ChEMBL for {len(ANTICANCER_TARGETS)} targets...")
    print(f"  (max {max_per_target} compounds per target, cached after first run)")
    data = fetch_real_data(
        targets=ANTICANCER_TARGETS,
        max_per_target=max_per_target,
        use_cache=use_cache,
    )

    print(f"\n  Enriching compounds with clinical phase data...")
    enrich_compounds_with_phase(data, use_cache=use_cache)

    save_data(data, DATA_PATH)
    return data


def merge_ground_truth_compounds(live_data: dict) -> dict:
    """
    Ensure all 18 ground-truth compounds are present in the dataset.

    The live ChEMBL fetch may not return all 18 curated compounds
    (e.g., if they have unusual assay types or are excluded by filters).
    We merge them in from the curated dataset to ensure the benchmark
    has all positive examples.
    """
    curated = load_curated_data()

    live_compound_ids = {c["chembl_id"] for c in live_data["compounds"]}
    live_pairs = {
        (b["compound"], b["target"])
        for b in live_data["bioactivities"]
    }

    added_compounds = 0
    added_activities = 0

    for compound in curated["compounds"]:
        cid = compound["chembl_id"]
        if cid.startswith("CHEMBL_"):
            continue  # skip test compounds
        if cid not in live_compound_ids:
            live_data["compounds"].append(compound)
            live_compound_ids.add(cid)
            added_compounds += 1

    for activity in curated["bioactivities"]:
        pair = (activity["compound"], activity["target"])
        if pair not in live_pairs:
            live_data["bioactivities"].append(activity)
            live_pairs.add(pair)
            added_activities += 1

    # Merge toxicity data (live data has empty toxicity)
    curated_tox = curated.get("toxicity", {})
    for cid, tox in curated_tox.items():
        if cid not in live_data.get("toxicity", {}):
            live_data.setdefault("toxicity", {})[cid] = tox

    if added_compounds or added_activities:
        print(f"  Merged {added_compounds} curated compounds, "
              f"{added_activities} curated bioactivities")

    return live_data


def run_scaleup_benchmark():
    """Run the full scale-up benchmark."""
    print("=" * 80)
    print(" SCALE-UP VALIDATION: ChemPath on Real ChEMBL Data")
    print("=" * 80)

    # Step 1: Get data
    print("\n[1] Data Acquisition")
    t0 = time.time()
    data = fetch_or_load_data(max_per_target=200)
    data = merge_ground_truth_compounds(data)
    t_data = time.time() - t0

    n_compounds = len(data["compounds"])
    n_targets = len(data["targets"])
    n_activities = len(data["bioactivities"])
    print(f"  Dataset: {n_compounds} compounds, {n_targets} targets, "
          f"{n_activities} bioactivities")
    print(f"  Data load time: {t_data:.1f}s")

    # Check ground truth coverage
    compound_ids = {c["chembl_id"] for c in data["compounds"]}
    gt_present = sum(1 for cid in GROUND_TRUTH if cid in compound_ids)
    print(f"  Ground-truth drugs in dataset: {gt_present}/{len(GROUND_TRUTH)}")

    # Step 2: Build graph
    print("\n[2] Building Multi-hop Graph")
    t0 = time.time()
    G = build_multihop_graph(data, verbose=True)
    t_graph = time.time() - t0
    summary = get_multihop_summary(G)
    print(f"  Build time: {t_graph:.1f}s")

    # Step 3: Baseline comparison
    print("\n[3] Benchmark: Recovering FDA-approved indications")
    k_values = [1, 3, 5, 10]

    print("  Running multi-hop benchmark...")
    t0 = time.time()
    multihop = run_multihop_benchmark(G, weight_key="w_efficacy", k_values=k_values)
    t_multihop = time.time() - t0

    print("  Running single-hop baseline...")
    single_hop = run_single_hop_baseline(G, k_values=k_values)

    print("  Running random baseline...")
    random_bl = run_random_baseline(G, k_values=k_values)

    report = format_benchmark_report(
        [multihop, single_hop, random_bl], show_per_drug=True
    )
    print(report)
    print(f"  Multi-hop pathfinding time: {t_multihop:.1f}s")

    # Step 4: AUROC / AUPR
    print("\n[4] AUROC / AUPR")
    t0 = time.time()
    metrics = compute_auroc_aupr(G, weight_key="w_efficacy")
    t_auroc = time.time() - t0

    print(f"  AUROC:          {metrics.auroc:.3f}  (random = 0.500)")
    print(f"  AUPR:           {metrics.aupr:.3f}  (random = {metrics.baseline_aupr:.3f})")
    print(f"  Positive pairs: {metrics.n_positive}")
    print(f"  Negative pairs: {metrics.n_negative}")
    print(f"  Total scored:   {metrics.n_total}")
    print(f"  Compute time:   {t_auroc:.1f}s")

    # Step 5: Comparison with curated baseline
    print("\n[5] Curated (18-drug) vs Scale-up Comparison")
    curated = load_curated_data()
    G_curated = build_multihop_graph(curated, verbose=False)
    m_curated = run_multihop_benchmark(G_curated, weight_key="w_efficacy", k_values=k_values)
    c_metrics = compute_auroc_aupr(G_curated, weight_key="w_efficacy")

    print(f"  {'Metric':<20s} | {'Curated (18)':>14s} | {'Scale-up ({0})':>14s}".format(
        n_compounds))
    print(f"  {'─' * 55}")
    print(f"  {'MRR':<20s} | {m_curated.mrr:>14.3f} | {multihop.mrr:>14.3f}")
    print(f"  {'Recall@1':<20s} | {m_curated.recall_at_k[1]:>13.1%} | {multihop.recall_at_k[1]:>13.1%}")
    print(f"  {'Recall@3':<20s} | {m_curated.recall_at_k[3]:>13.1%} | {multihop.recall_at_k[3]:>13.1%}")
    print(f"  {'AUROC':<20s} | {c_metrics.auroc:>14.3f} | {metrics.auroc:>14.3f}")
    print(f"  {'AUPR':<20s} | {c_metrics.aupr:>14.3f} | {metrics.aupr:>14.3f}")
    print(f"  {'Compounds':<20s} | {G_curated.number_of_nodes():>14d} | {G.number_of_nodes():>14d}")
    print(f"  {'Edges':<20s} | {G_curated.number_of_edges():>14d} | {G.number_of_edges():>14d}")

    # Verdict
    print(f"\n{'=' * 80}")
    print(" VERDICT")
    print(f"{'=' * 80}")
    if metrics.auroc >= 0.70:
        print(f"  AUROC = {metrics.auroc:.3f} >= 0.70  →  PASS: Framework scales.")
        print(f"  Multi-hop ranking maintains discriminative power with "
              f"{n_compounds} compounds.")
    elif metrics.auroc >= 0.60:
        print(f"  AUROC = {metrics.auroc:.3f}  →  MARGINAL: Moderate signal.")
        print(f"  Framework shows some scaling ability but may need "
              f"parameter tuning.")
    else:
        print(f"  AUROC = {metrics.auroc:.3f} < 0.60  →  FAIL: Signal lost at scale.")
        print(f"  Multi-hop ranking does not generalize beyond curated data.")

    if multihop.mrr > random_bl.mrr * 1.5:
        print(f"  MRR = {multihop.mrr:.3f} vs random {random_bl.mrr:.3f}  →  "
              f"{multihop.mrr / random_bl.mrr:.1f}x lift over random.")
    else:
        print(f"  MRR = {multihop.mrr:.3f} vs random {random_bl.mrr:.3f}  →  "
              f"Minimal lift over random.")
    print(f"{'=' * 80}")


if __name__ == "__main__":
    run_scaleup_benchmark()
