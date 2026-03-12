"""
Retrospective validation — tests whether ChemPath ranks known approved drugs highly.

This provides an honest self-assessment: given a target with multiple compounds,
do approved drugs (Phase 4) end up in the top rankings?

Key metric: Enrichment Factor — how much better than random is ChemPath at
finding approved drugs in the top N recommendations?
"""

import networkx as nx
from dataclasses import dataclass
from chempath.graph.optimizer import rank_compounds_for_target
from chempath.chemistry.properties import compute_enriched_properties


@dataclass
class ValidationResult:
    target_id: str
    target_name: str
    total_compounds: int
    approved_compounds: list[dict]
    top_n: int
    approved_in_top_n: list[dict]
    enrichment_factor: float
    hit_rate: float
    ranking_details: list[dict]


def validate_target(
    G: nx.DiGraph,
    target_id: str,
    top_n: int = 10,
    strategy: str = "balanced",
) -> ValidationResult | None:
    """Check if approved drugs rank in the top N for a specific target."""
    if target_id not in G:
        return None

    recs = rank_compounds_for_target(G, target_id, strategy=strategy)
    if not recs:
        return None

    total = len(recs)
    approved = []
    for rec in recs:
        phase = G.nodes.get(rec.compound_id, {}).get("phase", 0)
        if phase >= 4:
            approved.append({
                "name": rec.compound_name,
                "id": rec.compound_id,
                "rank": rec.rank,
                "ic50": rec.ic50_nm,
                "phase": phase,
            })

    if not approved:
        return None

    approved_in_top = [a for a in approved if a["rank"] <= top_n]
    n_approved = len(approved)
    expected_random = n_approved / total * top_n if total > 0 else 0
    ef = len(approved_in_top) / expected_random if expected_random > 0 else 0.0
    hit_rate = len(approved_in_top) / min(top_n, total) if total > 0 else 0.0

    ranking_details = []
    for rec in recs[:top_n]:
        phase = G.nodes.get(rec.compound_id, {}).get("phase", 0)
        ranking_details.append({
            "rank": rec.rank,
            "name": rec.compound_name,
            "ic50": rec.ic50_nm,
            "weight": rec.weight,
            "phase": phase,
            "is_approved": phase >= 4,
            "confidence": rec.confidence,
        })

    return ValidationResult(
        target_id=target_id,
        target_name=G.nodes[target_id].get("name", target_id),
        total_compounds=total,
        approved_compounds=approved,
        top_n=min(top_n, total),
        approved_in_top_n=approved_in_top,
        enrichment_factor=round(ef, 2),
        hit_rate=round(hit_rate, 3),
        ranking_details=ranking_details,
    )


def run_retrospective_validation(
    G: nx.DiGraph,
    top_n: int = 10,
    strategies: list[str] | None = None,
) -> dict[str, list[ValidationResult]]:
    """Run validation across all targets and strategies."""
    if strategies is None:
        strategies = ["balanced", "efficacy", "safety"]

    targets = [n for n, d in G.nodes(data=True) if d.get("node_type") == "target"]
    results = {}

    for strategy in strategies:
        strategy_results = []
        for tid in targets:
            result = validate_target(G, tid, top_n=top_n, strategy=strategy)
            if result is not None:
                strategy_results.append(result)
        results[strategy] = strategy_results

    return results


def compute_multi_objective_scores(
    G: nx.DiGraph,
    target_id: str,
) -> list[dict]:
    """
    Compute enriched multi-dimensional scores for all compounds targeting a protein.
    Returns data suitable for multi-objective Pareto analysis with >2 dimensions.
    """
    if target_id not in G:
        return []

    scores = []
    for compound_id in G.predecessors(target_id):
        node = G.nodes[compound_id]
        if node.get("node_type") != "compound":
            continue

        edge = G.edges[compound_id, target_id]
        smiles = node.get("smiles", "")
        if not smiles:
            continue

        props = compute_enriched_properties(smiles)

        scores.append({
            "compound_id": compound_id,
            "compound_name": node.get("name", compound_id),
            "phase": node.get("phase", 0),
            "ic50_nm": edge["ic50_nm"],
            "weight": edge["weight"],
            "source": edge.get("source", "experimental"),
            # Multi-dimensional properties
            "mw": props.estimated_mw,
            "logp": props.estimated_logp,
            "hbd": props.estimated_hbd,
            "hba": props.estimated_hba,
            "tpsa": props.estimated_tpsa,
            "rotatable_bonds": props.estimated_rotatable,
            "aromatic_rings": props.aromatic_rings,
            "qed": props.qed_score,
            "sa_score": props.sa_score,
            "lipinski_violations": props.lipinski_violations,
            "risk_score": props.risk_score,
        })

    return scores


def compute_multi_objective_pareto(scores: list[dict]) -> list[dict]:
    """
    Compute Pareto front over multiple objectives:
    Minimize: IC50, risk_score, sa_score (synthetic difficulty)
    Maximize: QED (convert to minimize by using 1-QED)

    A point is Pareto-optimal if no other point is strictly better on ALL objectives.
    """
    objectives = [
        ("ic50_nm", "minimize"),
        ("risk_score", "minimize"),
        ("sa_score", "minimize"),
        ("qed", "maximize"),
    ]

    def dominates(a, b):
        """Return True if a dominates b (a is at least as good on all, strictly better on one)."""
        at_least_as_good = True
        strictly_better = False
        for key, direction in objectives:
            va, vb = a[key], b[key]
            if direction == "maximize":
                va, vb = -va, -vb
            if va > vb:
                at_least_as_good = False
                break
            if va < vb:
                strictly_better = True
        return at_least_as_good and strictly_better

    pareto = []
    for candidate in scores:
        dominated = any(dominates(other, candidate) for other in scores if other is not candidate)
        if not dominated:
            pareto.append(candidate)

    pareto.sort(key=lambda x: x["ic50_nm"])
    return pareto


def format_validation_report(results: dict[str, list[ValidationResult]]) -> str:
    """Format validation results as a readable report."""
    lines = [
        f"\n{'=' * 70}",
        " RETROSPECTIVE VALIDATION: Do approved drugs rank highly?",
        f"{'=' * 70}",
    ]

    for strategy, vresults in results.items():
        lines.append(f"\n  Strategy: {strategy.upper()}")
        lines.append(f"  {'─' * 60}")

        if not vresults:
            lines.append("    No targets with approved drugs found.")
            continue

        total_ef = 0
        for vr in vresults:
            n_found = len(vr.approved_in_top_n)
            n_total = len(vr.approved_compounds)
            lines.append(
                f"    {vr.target_name:30s} | "
                f"Approved: {n_found}/{n_total} in top {vr.top_n} | "
                f"EF={vr.enrichment_factor:.1f}x | "
                f"Hit rate={vr.hit_rate:.1%}"
            )
            total_ef += vr.enrichment_factor

        avg_ef = total_ef / len(vresults) if vresults else 0
        lines.append(f"\n    Average Enrichment Factor: {avg_ef:.1f}x (>1.0 = better than random)")

    lines.append(f"\n{'=' * 70}")
    return "\n".join(lines)
