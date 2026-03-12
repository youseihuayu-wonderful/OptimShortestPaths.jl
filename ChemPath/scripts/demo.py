"""
ChemPath Demo — Full Pipeline
Runs: data → validate → graph → optimize → sensitivity analysis.
"""

from chempath.data.curated_data import load_curated_data
from chempath.chemistry.smiles import validate_batch
from chempath.graph.builder import build_drug_target_graph, get_graph_summary
from chempath.graph.optimizer import (
    rank_compounds_for_target, compute_pareto_front, find_knee_point,
    format_recommendations, format_confidence_detail,
)
from chempath.graph.analysis import run_full_sensitivity_analysis, format_sensitivity_report


def main():
    print("=" * 70)
    print(" ChemPath — Drug Screening Pipeline Demo")
    print("=" * 70)

    # Step 1: Load Data
    print("\n[Step 1] Loading curated ChEMBL data...")
    data = load_curated_data()
    print(f"  Loaded {len(data['compounds'])} compounds, {len(data['targets'])} targets, "
          f"{len(data['bioactivities'])} bioactivities")

    # Step 2: SMILES Validation
    print("\n[Step 2] Validating SMILES strings...")
    all_smiles = [c["smiles"] for c in data["compounds"]]
    batch_result = validate_batch(all_smiles)
    s = batch_result["summary"]
    print(f"  Results: {s['valid_count']}/{s['total']} valid ({s['valid_rate']:.0%})")
    if batch_result["invalid"]:
        print(f"  Rejected:")
        for inv in batch_result["invalid"]:
            print(f"    - '{inv.smiles[:40]}...' : {inv.error}")

    # Step 3: Build Graph
    print("\n[Step 3] Building drug-target graph...")
    G = build_drug_target_graph(data, toxicity_penalty=1.0)
    summary = get_graph_summary(G)
    print(f"  Nodes: {summary['total_nodes']} ({summary['compounds']} compounds, {summary['targets']} targets)")
    print(f"  Edges: {summary['total_edges']} ({summary['experimental_edges']} experimental)")

    # Step 4: Rank Compounds
    targets_to_analyze = [
        ("CHEMBL203", "EGFR"),
        ("CHEMBL5145", "BRAF"),
    ]

    for target_id, target_name in targets_to_analyze:
        print(f"\n[Step 4] Ranking compounds for {target_name}...")

        recs = rank_compounds_for_target(G, target_id, data.get("toxicity"), strategy="balanced")
        print(format_recommendations(recs, f"{target_name} Inhibitors — Balanced Strategy"))

        for rec in recs[:2]:
            print(format_confidence_detail(rec))

        # Step 5: Pareto Front
        print(f"\n[Step 5] Computing Pareto front for {target_name}...")
        pareto = compute_pareto_front(recs)
        print(format_recommendations(pareto, f"{target_name} — Pareto-Optimal Solutions"))

        knee = find_knee_point(pareto)
        if knee:
            print(f"\n  >>> KNEE POINT RECOMMENDATION: {knee.compound_name} <<<")
            print(f"      Best balance of efficacy (IC50={knee.ic50_nm}nM) and safety (Tox={knee.toxicity:.2f})")

    # Step 6: Sensitivity Analysis
    print(f"\n[Step 6] Running sensitivity analysis for EGFR...")
    results = run_full_sensitivity_analysis(data, "CHEMBL203")
    for result in results:
        print(format_sensitivity_report(result))

    # Summary
    print("\n" + "=" * 70)
    print(" PIPELINE SUMMARY")
    print("=" * 70)
    print(f"  Python:        3.12 (managed by uv)")
    print(f"  Data:          {s['valid_count']} valid compounds from {len(data['targets'])} targets")
    print(f"  Graph:         {summary['total_nodes']} nodes, {summary['total_edges']} edges")
    print(f"  Strategies:    Balanced | Efficacy-First | Safety-First")
    print(f"  Pareto:        Multi-objective optimization (efficacy vs toxicity)")
    print(f"  Confidence:    HIGH / MEDIUM / LOW labels on every recommendation")
    print(f"  Sensitivity:   3 analyses (toxicity weight, IC50 noise, missing data)")
    print("=" * 70)


if __name__ == "__main__":
    main()
