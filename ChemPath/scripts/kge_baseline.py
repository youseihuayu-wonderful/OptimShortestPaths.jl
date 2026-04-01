"""
KGE Baseline — TransE + ComplEx on Hetionet for drug repurposing AUROC.

Uses PyKEEN to train TransE and ComplEx on Hetionet edges (excluding
'treats' edges which are held out for evaluation). Evaluates on
PharmacotherapyDB indications using the same protocol as our
shortest-path method.

ALL DATA IS REAL:
  - Training triples: Hetionet v1.0 edges (minus 'treats')
  - Evaluation: PharmacotherapyDB 755 disease-modifying indications
  - Same compound/disease sets as our benchmark
"""
from __future__ import annotations

import json
import math
import os
import pickle
import sys
import time
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))

from chempath_enriched_benchmark import (
    load_hetionet_enriched,
    load_pharmacotherapydb,
)

CACHE_DIR = Path(__file__).resolve().parent.parent / "data" / "hetionet"
RESULTS_DIR = Path(__file__).resolve().parent.parent / "data" / "paper_results"


def build_triples(G, ground_truth):
    """Build (head, relation, tail) triples from Hetionet graph.

    Excludes 'treats' edges (held out for evaluation).
    Uses REAL graph edges only.
    """
    triples = []
    for u, v, data in G.edges(data=True):
        edge_type = data.get("edge_type", "unknown")
        triples.append((u, edge_type, v))
    # ground_truth contains the treats edges (held out during graph construction)
    # These are NOT included in training — they're what we evaluate on
    print(f"  Training triples: {len(triples)} (excludes {sum(len(v) for v in ground_truth.values())} treats edges)")
    return triples


def run_pykeen_model(triples, model_name, eval_pairs, positive_pairs, epochs=100):
    """Train a PyKEEN model and evaluate on drug repurposing.

    Returns (auroc, auprc, training_time).
    """
    import numpy as np
    import torch
    from pykeen.triples import TriplesFactory
    from pykeen.pipeline import pipeline

    print(f"\n  Training {model_name} ({epochs} epochs)...")

    # Build triples factory — PyKEEN requires numpy array of shape (n, 3)
    triple_array = np.array([[h, r, t] for h, r, t in triples], dtype=str)
    tf = TriplesFactory.from_labeled_triples(
        triples=triple_array,
    )

    # Split into train/test (90/10) for PyKEEN pipeline requirement
    # We don't use PyKEEN's test evaluation — we evaluate on PharmacotherapyDB
    train_tf, test_tf = tf.split(ratios=[0.9, 0.1], random_state=42)
    print(f"  Train triples: {train_tf.num_triples}, Test triples: {test_tf.num_triples}")

    # Train
    t0 = time.time()
    result = pipeline(
        training=train_tf,
        testing=test_tf,
        model=model_name,
        training_kwargs=dict(
            num_epochs=epochs,
            batch_size=1024,
        ),
        model_kwargs=dict(
            embedding_dim=128,
        ),
        optimizer_kwargs=dict(
            lr=0.001,
        ),
        random_seed=42,
        device="cpu",
    )
    train_time = time.time() - t0
    print(f"  Trained in {train_time:.1f}s")

    model = result.model
    model.eval()

    # Score all compound-disease pairs
    print(f"  Scoring {len(eval_pairs)} compound-disease pairs...")
    scored_list = []

    for compound, disease in eval_pairs:
        label = 1 if (compound, disease) in positive_pairs else 0

        # PyKEEN score: model predicts (compound, treats, disease)
        try:
            h_id = tf.entity_to_id.get(compound)
            t_id = tf.entity_to_id.get(disease)
            # Use 'palliates' as closest relation to 'treats' (treats was held out)
            r_id = tf.relation_to_id.get("palliates")
            if r_id is None:
                # Fallback: use 'associates'
                r_id = tf.relation_to_id.get("associates")

            if h_id is not None and t_id is not None and r_id is not None:
                h_tensor = torch.tensor([[h_id]], dtype=torch.long)
                r_tensor = torch.tensor([[r_id]], dtype=torch.long)
                t_tensor = torch.tensor([[t_id]], dtype=torch.long)
                with torch.no_grad():
                    score = model.score_hrt(
                        torch.cat([h_tensor, r_tensor, t_tensor], dim=1)
                    ).item()
            else:
                score = -1e6  # Unknown entity
        except Exception:
            score = -1e6

        scored_list.append((score, label))

    # Compute AUROC
    scored_list.sort(key=lambda x: -x[0])
    n_pos = sum(lab for _, lab in scored_list)
    n_neg = len(scored_list) - n_pos

    if n_pos == 0 or n_neg == 0:
        return 0.5, 0.0, train_time

    rank_sum = 0.0
    for i, (score, label) in enumerate(scored_list):
        if label == 1:
            rank_sum += (len(scored_list) - i)
    auroc = (rank_sum - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)

    # AUPRC
    tp = 0
    auprc_sum = 0.0
    for i, (score, label) in enumerate(scored_list):
        if label == 1:
            tp += 1
            precision = tp / (i + 1)
            auprc_sum += precision
    auprc = auprc_sum / n_pos if n_pos > 0 else 0.0

    return auroc, auprc, train_time


def main():
    print("=" * 80)
    print("  KGE BASELINES — TransE + ComplEx on Hetionet")
    print("  ALL REAL DATA. No mock triples.")
    print("=" * 80)
    t_start = time.time()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Load data
    print("\n[1/3] Loading Hetionet...")
    het_data = load_hetionet_enriched()
    G = het_data["G"]
    ground_truth = het_data["ground_truth"]
    print(f"  Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    print("\n[2/3] Building training triples...")
    triples = build_triples(G, ground_truth)

    print("\n[3/3] Loading evaluation set...")
    positive_pairs, eval_compounds, eval_diseases = load_pharmacotherapydb(G, ground_truth)
    eval_pairs = [(c, d) for c in eval_compounds for d in eval_diseases]
    print(f"  Eval pairs: {len(eval_pairs)} ({len(eval_compounds)} compounds × {len(eval_diseases)} diseases)")
    print(f"  Positive: {len(positive_pairs)}")

    # Run models
    results = {}
    for model_name in ["TransE", "ComplEx"]:
        auroc, auprc, train_time = run_pykeen_model(
            triples, model_name, eval_pairs, positive_pairs, epochs=100
        )
        results[model_name] = {
            "auroc": auroc, "auprc": auprc, "train_time": train_time,
        }
        print(f"\n  {model_name}: AUROC={auroc:.4f}, AUPRC={auprc:.4f}, Time={train_time:.1f}s")

    # Summary comparison
    print("\n" + "=" * 80)
    print("  KGE BASELINE COMPARISON")
    print("=" * 80)
    print(f"  {'Method':<20s} {'AUROC':>8s} {'AUPRC':>8s} {'Time':>8s}")
    print(f"  {'─'*20} {'─'*8} {'─'*8} {'─'*8}")
    for m, r in results.items():
        print(f"  {m:<20s} {r['auroc']:>8.4f} {r['auprc']:>8.4f} {r['train_time']:>7.1f}s")

    # Save
    output_path = RESULTS_DIR / "kge_baseline.json"
    with open(output_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved to: {output_path}")

    print(f"\n  Total time: {time.time() - t_start:.0f}s")


if __name__ == "__main__":
    main()
