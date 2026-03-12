"""
Retrospective benchmark: Can ChemPath recover known drug-disease pairs?

Ground truth: FDA-approved drug indications from DailyMed/DrugBank.
Test: For each approved drug, run multi-hop pathfinding and check whether its
real-world indication(s) appear in the top-ranked diseases.

Metrics:
  - Recall@k: fraction of true indications found in top k results
  - Mean Reciprocal Rank (MRR): 1/rank of the first true indication
  - Enrichment Factor (EF): how much better than random
  - AUROC: area under ROC curve (all compound-disease pairs)
  - AUPR: area under precision-recall curve

Baselines:
  - Random: diseases ranked randomly
  - Single-hop IC50: rank diseases by best IC50 of any compound→target edge
    (ignores pathway/disease layers — just sorts by binding affinity)
"""

import math
import random
import networkx as nx
from dataclasses import dataclass, field
from collections import defaultdict

from chempath.graph.pathfinding import find_repurposing_candidates, scan_all_repurposing


# ---------------------------------------------------------------------------
# Ground truth: Known FDA-approved drug → disease associations
# ---------------------------------------------------------------------------

# These are real-world approved indications from FDA labels.
# compound_id → list of disease_ids in our network.
GROUND_TRUTH = {
    # --- EGFR inhibitors → NSCLC ---
    "CHEMBL939": ["DIS_NSCLC"],        # Gefitinib
    "CHEMBL553": ["DIS_NSCLC"],        # Erlotinib
    "CHEMBL1201585": ["DIS_NSCLC"],    # Osimertinib
    # --- Dual EGFR/HER2 → Breast Cancer ---
    "CHEMBL554": ["DIS_BREAST"],       # Lapatinib
    "CHEMBL180022": ["DIS_BREAST"],    # Neratinib
    # --- ALK inhibitors → NSCLC (ALK-rearranged) ---
    "CHEMBL3545252": ["DIS_NSCLC"],    # Alectinib
    "CHEMBL2007641": ["DIS_NSCLC"],    # Ceritinib
    "CHEMBL601719": ["DIS_NSCLC"],     # Crizotinib (ALK + MET)
    # --- BRAF inhibitors → Melanoma ---
    "CHEMBL1229517": ["DIS_MELANOMA"], # Vemurafenib
    "CHEMBL2028663": ["DIS_MELANOMA"], # Dabrafenib
    # --- BCR-ABL inhibitors → CML ---
    "CHEMBL941": ["DIS_CML", "DIS_GIST"],  # Imatinib (ABL1 + KIT)
    "CHEMBL288441": ["DIS_CML"],       # Bosutinib
    "CHEMBL255863": ["DIS_CML"],       # Nilotinib
    "CHEMBL1171837": ["DIS_CML"],      # Ponatinib
    "CHEMBL1421": ["DIS_CML"],         # Dasatinib
    # --- Multi-kinase ---
    "CHEMBL1336": ["DIS_HCC", "DIS_RCC", "DIS_THYROID"],  # Sorafenib (VEGFR2 + BRAF)
    "CHEMBL24828": ["DIS_THYROID"],    # Vandetanib (VEGFR2 + EGFR)
    # --- mTOR ---
    "CHEMBL1642": ["DIS_RCC", "DIS_BREAST"],  # Everolimus
}


@dataclass
class BenchmarkResult:
    """Result for a single drug in the benchmark."""
    compound_id: str
    compound_name: str
    true_diseases: list[str]
    predicted_diseases: list[str]  # ranked list of disease IDs
    predicted_probs: list[float]   # probabilities for each prediction
    rank_of_first_hit: int | None  # rank of first true disease (1-indexed), None if not found
    hits_in_top_k: dict[int, int]  # {k: n_hits_in_top_k}
    all_ranks: list[int]           # ranks of all true diseases found


@dataclass
class BenchmarkSummary:
    """Aggregate metrics across all benchmark drugs."""
    method_name: str
    n_drugs: int
    n_drugs_with_hits: int
    recall_at_k: dict[int, float]     # {k: recall}
    mrr: float                        # mean reciprocal rank
    enrichment_factor: float          # vs random baseline
    per_drug_results: list[BenchmarkResult]
    n_total_diseases: int


# ---------------------------------------------------------------------------
# Benchmark runners
# ---------------------------------------------------------------------------

def run_multihop_benchmark(
    G: nx.DiGraph,
    weight_key: str = "w_efficacy",
    k_values: list[int] | None = None,
) -> BenchmarkSummary:
    """
    Run the multi-hop pathfinding benchmark.

    For each ground-truth drug, find repurposing candidates and check
    whether the known indication(s) appear in the ranked results.
    """
    if k_values is None:
        k_values = [1, 3, 5, 10]

    disease_nodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "disease"]
    n_diseases = len(disease_nodes)

    results = []
    for compound_id, true_diseases in GROUND_TRUTH.items():
        if compound_id not in G:
            continue

        compound_name = G.nodes[compound_id].get("name", compound_id)

        # Run repurposing search
        repurposing = find_repurposing_candidates(
            G, compound_id, k_paths=3, weight_key=weight_key
        )

        predicted_diseases = [r.disease_id for r in repurposing]
        predicted_probs = [r.best_probability for r in repurposing]

        # Find ranks of true diseases
        all_ranks = []
        for true_dis in true_diseases:
            if true_dis in predicted_diseases:
                rank = predicted_diseases.index(true_dis) + 1
                all_ranks.append(rank)

        rank_first = min(all_ranks) if all_ranks else None

        hits_in_top_k = {}
        for k in k_values:
            top_k_set = set(predicted_diseases[:k])
            hits = len(set(true_diseases) & top_k_set)
            hits_in_top_k[k] = hits

        results.append(BenchmarkResult(
            compound_id=compound_id,
            compound_name=compound_name,
            true_diseases=true_diseases,
            predicted_diseases=predicted_diseases,
            predicted_probs=predicted_probs,
            rank_of_first_hit=rank_first,
            hits_in_top_k=hits_in_top_k,
            all_ranks=all_ranks,
        ))

    return _compute_summary("Multi-hop pathfinding", results, k_values, n_diseases)


def run_random_baseline(
    G: nx.DiGraph,
    n_trials: int = 100,
    k_values: list[int] | None = None,
    seed: int = 42,
) -> BenchmarkSummary:
    """
    Random baseline: diseases ranked randomly.
    Averaged over n_trials to get stable estimates.
    """
    if k_values is None:
        k_values = [1, 3, 5, 10]

    disease_nodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "disease"]
    n_diseases = len(disease_nodes)
    rng = random.Random(seed)

    # Average metrics over trials
    all_mrrs = []
    all_recall_at_k = {k: [] for k in k_values}

    for _ in range(n_trials):
        trial_results = []
        for compound_id, true_diseases in GROUND_TRUTH.items():
            if compound_id not in G:
                continue

            shuffled = disease_nodes.copy()
            rng.shuffle(shuffled)

            all_ranks = []
            for true_dis in true_diseases:
                if true_dis in shuffled:
                    rank = shuffled.index(true_dis) + 1
                    all_ranks.append(rank)

            rank_first = min(all_ranks) if all_ranks else None

            hits_in_top_k = {}
            for k in k_values:
                top_k_set = set(shuffled[:k])
                hits = len(set(true_diseases) & top_k_set)
                hits_in_top_k[k] = hits

            trial_results.append(BenchmarkResult(
                compound_id=compound_id,
                compound_name=G.nodes[compound_id].get("name", compound_id),
                true_diseases=true_diseases,
                predicted_diseases=shuffled,
                predicted_probs=[1.0 / n_diseases] * n_diseases,
                rank_of_first_hit=rank_first,
                hits_in_top_k=hits_in_top_k,
                all_ranks=all_ranks,
            ))

        # Compute trial metrics
        rrs = []
        for r in trial_results:
            if r.rank_of_first_hit:
                rrs.append(1.0 / r.rank_of_first_hit)
            else:
                rrs.append(0.0)
        all_mrrs.append(sum(rrs) / len(rrs) if rrs else 0.0)

        for k in k_values:
            total_true = sum(len(r.true_diseases) for r in trial_results)
            total_hits = sum(r.hits_in_top_k[k] for r in trial_results)
            all_recall_at_k[k].append(total_hits / total_true if total_true else 0.0)

    # Average across trials - build a representative result
    avg_mrr = sum(all_mrrs) / len(all_mrrs)
    avg_recall = {k: sum(v) / len(v) for k, v in all_recall_at_k.items()}

    # Return summary with averaged metrics
    summary = BenchmarkSummary(
        method_name="Random baseline",
        n_drugs=len([c for c in GROUND_TRUTH if c in G]),
        n_drugs_with_hits=0,  # not meaningful for averaged random
        recall_at_k=avg_recall,
        mrr=avg_mrr,
        enrichment_factor=1.0,  # by definition
        per_drug_results=[],
        n_total_diseases=n_diseases,
    )
    return summary


def run_single_hop_baseline(
    G: nx.DiGraph,
    k_values: list[int] | None = None,
) -> BenchmarkSummary:
    """
    Single-hop baseline: rank diseases by best IC50 of any compound→target
    edge, ignoring pathway and disease layers entirely.

    For each compound, find its targets, then rank diseases by the number
    of pathways shared (a proxy for single-hop reasoning without
    multi-hop path probabilities).
    """
    if k_values is None:
        k_values = [1, 3, 5, 10]

    disease_nodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "disease"]
    n_diseases = len(disease_nodes)

    results = []
    for compound_id, true_diseases in GROUND_TRUTH.items():
        if compound_id not in G:
            continue

        compound_name = G.nodes[compound_id].get("name", compound_id)

        # Find targets this compound hits
        targets = [
            n for n in G.successors(compound_id)
            if G.nodes[n].get("node_type") == "target"
        ]

        # Count diseases reachable via each target (unweighted, just counting edges)
        disease_scores = defaultdict(float)
        for tid in targets:
            ic50 = G.edges[compound_id, tid].get("ic50_nm", 1000)
            # Simple heuristic: score = 1/IC50 (higher is better)
            target_score = 1.0 / max(ic50, 0.01)

            # For each pathway the target reaches
            for pathway in G.successors(tid):
                if G.nodes[pathway].get("node_type") != "pathway":
                    continue
                for disease in G.successors(pathway):
                    if G.nodes[disease].get("node_type") != "disease":
                        continue
                    # Accumulate score (more targets/pathways = higher score)
                    disease_scores[disease] += target_score

        # Rank diseases by accumulated score
        ranked = sorted(disease_scores.keys(), key=lambda d: -disease_scores[d])
        # Add diseases not reached (at the end)
        for d in disease_nodes:
            if d not in ranked:
                ranked.append(d)

        all_ranks = []
        for true_dis in true_diseases:
            if true_dis in ranked:
                all_ranks.append(ranked.index(true_dis) + 1)

        rank_first = min(all_ranks) if all_ranks else None
        hits_in_top_k = {}
        for k in k_values:
            top_k_set = set(ranked[:k])
            hits = len(set(true_diseases) & top_k_set)
            hits_in_top_k[k] = hits

        results.append(BenchmarkResult(
            compound_id=compound_id,
            compound_name=compound_name,
            true_diseases=true_diseases,
            predicted_diseases=ranked,
            predicted_probs=[disease_scores.get(d, 0) for d in ranked],
            rank_of_first_hit=rank_first,
            hits_in_top_k=hits_in_top_k,
            all_ranks=all_ranks,
        ))

    return _compute_summary("Single-hop IC50 baseline", results, k_values, n_diseases)


# ---------------------------------------------------------------------------
# Metric computation
# ---------------------------------------------------------------------------

def _compute_summary(
    method_name: str,
    results: list[BenchmarkResult],
    k_values: list[int],
    n_diseases: int,
) -> BenchmarkSummary:
    """Compute aggregate metrics from per-drug results."""
    if not results:
        return BenchmarkSummary(
            method_name=method_name, n_drugs=0, n_drugs_with_hits=0,
            recall_at_k={k: 0.0 for k in k_values}, mrr=0.0,
            enrichment_factor=0.0, per_drug_results=[], n_total_diseases=n_diseases,
        )

    n_drugs = len(results)
    n_with_hits = sum(1 for r in results if r.rank_of_first_hit is not None)

    # Recall@k: fraction of true indications found in top k
    recall_at_k = {}
    for k in k_values:
        total_true = sum(len(r.true_diseases) for r in results)
        total_found = sum(r.hits_in_top_k.get(k, 0) for r in results)
        recall_at_k[k] = total_found / total_true if total_true > 0 else 0.0

    # MRR: mean reciprocal rank
    rrs = []
    for r in results:
        if r.rank_of_first_hit is not None:
            rrs.append(1.0 / r.rank_of_first_hit)
        else:
            rrs.append(0.0)
    mrr = sum(rrs) / len(rrs) if rrs else 0.0

    # Enrichment factor at k=3: how much better than random?
    # Random expected recall@3 = 3/n_diseases
    k_ef = 3
    random_recall = min(k_ef, n_diseases) / n_diseases if n_diseases > 0 else 0
    ef = recall_at_k.get(k_ef, 0) / random_recall if random_recall > 0 else 0.0

    return BenchmarkSummary(
        method_name=method_name,
        n_drugs=n_drugs,
        n_drugs_with_hits=n_with_hits,
        recall_at_k=recall_at_k,
        mrr=mrr,
        enrichment_factor=ef,
        per_drug_results=results,
        n_total_diseases=n_diseases,
    )


# ---------------------------------------------------------------------------
# Dimension comparison: efficacy vs safety vs evidence
# ---------------------------------------------------------------------------

def run_dimension_comparison(
    G: nx.DiGraph,
    k_values: list[int] | None = None,
) -> dict[str, BenchmarkSummary]:
    """
    Run the benchmark with each probability dimension as the weight key.
    Shows how results differ when optimizing for efficacy vs safety vs evidence.
    """
    dimensions = {
        "efficacy": "w_efficacy",
        "safety": "w_safety",
        "evidence": "w_evidence",
    }
    results = {}
    for name, weight_key in dimensions.items():
        results[name] = run_multihop_benchmark(G, weight_key=weight_key, k_values=k_values)
        results[name].method_name = f"Multi-hop ({name})"
    return results


# ---------------------------------------------------------------------------
# Report formatting
# ---------------------------------------------------------------------------

def format_benchmark_report(
    summaries: list[BenchmarkSummary],
    show_per_drug: bool = True,
) -> str:
    """Format benchmark results as a readable report."""
    lines = [
        f"\n{'=' * 80}",
        " RETROSPECTIVE BENCHMARK: Recovery of known drug-disease pairs",
        f"{'=' * 80}",
    ]

    # Summary table header
    k_values = sorted(summaries[0].recall_at_k.keys()) if summaries else []
    header = f"  {'Method':<30s} | {'MRR':>6s} | {'EF@3':>6s}"
    for k in k_values:
        header += f" | {'R@'+str(k):>6s}"
    lines.append(header)
    lines.append(f"  {'─' * (len(header) - 2)}")

    for s in summaries:
        row = f"  {s.method_name:<30s} | {s.mrr:>6.3f} | {s.enrichment_factor:>6.1f}"
        for k in k_values:
            row += f" | {s.recall_at_k[k]:>6.1%}"
        lines.append(row)

    lines.append(f"  {'─' * (len(header) - 2)}")
    lines.append(f"  Drugs evaluated: {summaries[0].n_drugs if summaries else 0}")
    lines.append(f"  Diseases in network: {summaries[0].n_total_diseases if summaries else 0}")
    lines.append(f"  EF@3 = Enrichment Factor at top 3 (>1.0 = better than random)")
    lines.append(f"  MRR = Mean Reciprocal Rank (1.0 = always rank 1)")

    if show_per_drug and summaries:
        # Show per-drug detail for the first (main) method
        main = summaries[0]
        if main.per_drug_results:
            lines.append(f"\n  Per-drug detail ({main.method_name}):")
            lines.append(f"  {'─' * 70}")
            for r in sorted(main.per_drug_results, key=lambda x: x.rank_of_first_hit or 999):
                true_str = ", ".join(r.true_diseases)
                if r.rank_of_first_hit is not None:
                    rank_str = f"rank={r.rank_of_first_hit}"
                    found = "HIT" if r.rank_of_first_hit <= 3 else "   "
                else:
                    rank_str = "NOT FOUND"
                    found = "MISS"
                lines.append(
                    f"    [{found}] {r.compound_name:25s} → {true_str:30s} | {rank_str}"
                )

    lines.append(f"\n{'=' * 80}")
    return "\n".join(lines)


def run_full_benchmark(G: nx.DiGraph, verbose: bool = True) -> dict:
    """
    Run the complete benchmark suite and return all results.

    Returns dict with keys: 'multihop', 'random', 'single_hop', 'dimensions', 'report'
    """
    k_values = [1, 3, 5, 10]

    # Main methods
    multihop = run_multihop_benchmark(G, weight_key="w_efficacy", k_values=k_values)
    random_bl = run_random_baseline(G, k_values=k_values)
    single_hop = run_single_hop_baseline(G, k_values=k_values)

    # Dimension comparison
    dimensions = run_dimension_comparison(G, k_values=k_values)

    # Format report
    report = format_benchmark_report(
        [multihop, single_hop, random_bl] + list(dimensions.values()),
        show_per_drug=True,
    )

    if verbose:
        print(report)

    return {
        "multihop": multihop,
        "random": random_bl,
        "single_hop": single_hop,
        "dimensions": dimensions,
        "report": report,
    }


# ---------------------------------------------------------------------------
# AUROC / AUPR: Binary classification metrics for drug-disease prediction
# ---------------------------------------------------------------------------

@dataclass
class ClassificationMetrics:
    """AUROC and AUPR results for binary drug-disease prediction."""
    auroc: float
    aupr: float
    n_positive: int
    n_negative: int
    n_total: int
    baseline_aupr: float  # random baseline = n_positive / n_total


def compute_auroc_aupr(
    G: nx.DiGraph,
    weight_key: str = "w_efficacy",
) -> ClassificationMetrics:
    """
    Compute AUROC and AUPR for drug-disease prediction.

    Treats every (compound, disease) pair as a binary classification problem:
      - Positive: pair is in GROUND_TRUTH
      - Negative: pair is NOT in GROUND_TRUTH
      - Score: path probability from multi-hop pathfinding (0 if unreachable)

    Returns ClassificationMetrics with AUROC, AUPR, and counts.
    """
    compounds = [n for n, d in G.nodes(data=True)
                 if d.get("node_type") == "compound"
                 and not n.startswith("CHEMBL_")]  # exclude test compounds
    diseases = [n for n, d in G.nodes(data=True)
                if d.get("node_type") == "disease"]

    # Build positive set
    positive_pairs = set()
    for cid, dids in GROUND_TRUTH.items():
        if cid in G:
            for did in dids:
                if did in G:
                    positive_pairs.add((cid, did))

    # Score all (compound, disease) pairs
    scores = []
    labels = []

    for cid in compounds:
        repurposing = find_repurposing_candidates(G, cid, k_paths=1, weight_key=weight_key)
        disease_probs = {r.disease_id: r.best_probability for r in repurposing}

        for did in diseases:
            score = disease_probs.get(did, 0.0)
            label = 1 if (cid, did) in positive_pairs else 0
            scores.append(score)
            labels.append(label)

    n_pos = sum(labels)
    n_neg = len(labels) - n_pos

    auroc = _compute_auroc(labels, scores)
    aupr = _compute_aupr(labels, scores)
    baseline_aupr = n_pos / len(labels) if labels else 0.0

    return ClassificationMetrics(
        auroc=auroc,
        aupr=aupr,
        n_positive=n_pos,
        n_negative=n_neg,
        n_total=len(labels),
        baseline_aupr=baseline_aupr,
    )


def _compute_auroc(labels: list[int], scores: list[float]) -> float:
    """
    Compute Area Under ROC Curve using the Wilcoxon-Mann-Whitney statistic.

    AUROC = P(score(positive) > score(negative)) over all pos-neg pairs.
    """
    n_pos = sum(labels)
    n_neg = len(labels) - n_pos
    if n_pos == 0 or n_neg == 0:
        return 0.0

    # Sort by score descending
    pairs = sorted(zip(scores, labels), key=lambda x: -x[0])

    # Count concordant pairs
    tp = 0
    concordant = 0.0
    for score, label in pairs:
        if label == 1:
            tp += 1
        else:
            # This negative is ranked below tp positives
            concordant += tp

    return concordant / (n_pos * n_neg)


def _compute_aupr(labels: list[int], scores: list[float]) -> float:
    """
    Compute Area Under Precision-Recall Curve.

    Uses the trapezoidal rule on the precision-recall curve.
    """
    n_pos = sum(labels)
    if n_pos == 0:
        return 0.0

    # Sort by score descending
    pairs = sorted(zip(scores, labels), key=lambda x: -x[0])

    tp = 0
    fp = 0
    aupr = 0.0
    prev_recall = 0.0

    for score, label in pairs:
        if label == 1:
            tp += 1
            precision = tp / (tp + fp)
            recall = tp / n_pos
            # Area of rectangle from prev_recall to recall at this precision
            aupr += precision * (recall - prev_recall)
            prev_recall = recall
        else:
            fp += 1

    return aupr
