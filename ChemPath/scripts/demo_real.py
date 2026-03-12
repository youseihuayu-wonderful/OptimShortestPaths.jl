"""
ChemPath Demo — Real ChEMBL Data Pipeline
Runs the full pipeline on real IC50 data fetched from ChEMBL.
"""

from pathlib import Path

from chempath.data.chembl_client import load_saved_data, ANTICANCER_TARGETS
from chempath.chemistry.smiles import validate_batch
from chempath.graph.builder import build_drug_target_graph, get_graph_summary
from chempath.graph.optimizer import (
    rank_compounds_for_target, compute_pareto_front, find_knee_point,
    format_recommendations, format_confidence_detail,
)

DATA_PATH = Path(__file__).parent.parent / "data" / "chembl_real.json"


def main():
    print("=" * 70)
    print(" ChemPath — Real ChEMBL Data Pipeline")
    print("=" * 70)

    # Step 1: Load saved data
    print("\n[Step 1] Loading real ChEMBL data...")
    data = load_saved_data(DATA_PATH)
    print(f"  {len(data['compounds'])} compounds, {len(data['targets'])} targets, "
          f"{len(data['bioactivities'])} bioactivities")

    # Step 2: SMILES Validation
    print("\n[Step 2] Validating SMILES strings...")
    all_smiles = [c["smiles"] for c in data["compounds"]]
    batch_result = validate_batch(all_smiles)
    s = batch_result["summary"]
    print(f"  Results: {s['valid_count']}/{s['total']} valid ({s['valid_rate']:.1%})")
    if batch_result["invalid"]:
        print(f"  Rejected {s['invalid_count']} compounds with invalid SMILES")

    # Step 3: Build Graph
    print("\n[Step 3] Building drug-target graph...")
    G = build_drug_target_graph(data, toxicity_penalty=0.0, verbose=True)
    summary = get_graph_summary(G)
    print(f"  Nodes: {summary['total_nodes']} "
          f"({summary['compounds']} compounds, {summary['targets']} targets)")
    print(f"  Edges: {summary['total_edges']} ({summary['experimental_edges']} experimental)")

    # Step 4: Rank compounds for each target
    targets_to_analyze = [
        ("CHEMBL203", "EGFR"),
        ("CHEMBL5145", "BRAF"),
        ("CHEMBL4282", "ALK"),
        ("CHEMBL2842", "VEGFR2"),
    ]

    for target_id, target_name in targets_to_analyze:
        print(f"\n[Step 4] Ranking compounds for {target_name}...")
        recs = rank_compounds_for_target(G, target_id, strategy="balanced")

        if not recs:
            print(f"  No compounds found for {target_name}")
            continue

        # Show top 5
        top5 = recs[:5]
        print(format_recommendations(top5, f"Top 5 {target_name} Inhibitors — Balanced"))

        for rec in top5[:2]:
            print(format_confidence_detail(rec))

        # Step 5: Pareto Front
        print(f"\n[Step 5] Pareto front for {target_name}...")
        pareto = compute_pareto_front(recs)
        print(f"  {len(pareto)} Pareto-optimal solutions from {len(recs)} candidates")
        if pareto:
            top_pareto = pareto[:5]
            print(format_recommendations(top_pareto, f"{target_name} — Pareto-Optimal (top 5)"))

            knee = find_knee_point(pareto)
            if knee:
                print(f"\n  >>> KNEE POINT: {knee.compound_name} "
                      f"(IC50={knee.ic50_nm:.1f}nM, Weight={knee.weight:.4f}) <<<")

    # Summary
    print("\n" + "=" * 70)
    print(" PIPELINE SUMMARY — REAL DATA")
    print("=" * 70)
    print(f"  Source:        ChEMBL REST API")
    print(f"  Compounds:     {summary['compounds']} (SMILES-validated)")
    print(f"  Targets:       {summary['targets']}")
    print(f"  Edges:         {summary['total_edges']}")
    print(f"  Targets analyzed: {len(targets_to_analyze)}")
    print("=" * 70)


if __name__ == "__main__":
    main()
