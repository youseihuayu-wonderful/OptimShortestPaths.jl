"""
Multi-objective drug screening with shortest-path optimization and confidence labels.

Features:
  - Single-target compound ranking (3 strategies)
  - Multi-objective Pareto front (efficacy vs toxicity)
  - Confidence labeling: HIGH / MEDIUM / LOW per recommendation
"""

import math
import networkx as nx
from dataclasses import dataclass, field
from chempath.chemistry.properties import estimate_qed, estimate_synthetic_accessibility


@dataclass
class DrugRecommendation:
    """A single drug recommendation with confidence label."""
    rank: int
    compound_id: str
    compound_name: str
    target_id: str
    target_name: str
    ic50_nm: float
    weight: float
    toxicity: float
    source: str
    confidence: str
    confidence_reasons: list[str] = field(default_factory=list)


def assign_confidence(
    ic50_nm: float, source: str, phase: int, toxicity: float
) -> tuple[str, list[str]]:
    """
    Assign confidence label based on 4 dimensions:
    data source, clinical phase, potency, toxicity.

    HIGH (score>=8): experimental + approved + potent + low-tox
    MEDIUM (5-7): experimental but early phase or moderate tox
    LOW (<5): predicted, or high tox, or no clinical data
    """
    reasons = []
    score = 0

    if source == "experimental":
        score += 3
        reasons.append("Experimental IC50 data")
    else:
        reasons.append("AI-predicted IC50 (not experimentally verified)")

    if phase >= 4:
        score += 3
        reasons.append(f"Approved drug (Phase {phase})")
    elif phase >= 3:
        score += 2
        reasons.append(f"Phase {phase} clinical trial")
    elif phase >= 1:
        score += 1
        reasons.append(f"Phase {phase} clinical trial")
    else:
        reasons.append("No clinical trial data")

    if ic50_nm < 10:
        score += 2
        reasons.append(f"High potency (IC50={ic50_nm}nM)")
    elif ic50_nm < 100:
        score += 1
        reasons.append(f"Moderate potency (IC50={ic50_nm}nM)")
    else:
        reasons.append(f"Low potency (IC50={ic50_nm}nM)")

    if toxicity < 0.2:
        score += 2
        reasons.append(f"Low toxicity ({toxicity:.2f})")
    elif toxicity < 0.35:
        score += 1
        reasons.append(f"Moderate toxicity ({toxicity:.2f})")
    else:
        reasons.append(f"High toxicity ({toxicity:.2f})")

    if score >= 8:
        return "HIGH", reasons
    elif score >= 5:
        return "MEDIUM", reasons
    else:
        return "LOW", reasons


def rank_compounds_for_target(
    G: nx.DiGraph,
    target_id: str,
    toxicity_data: dict | None = None,
    strategy: str = "balanced",
) -> list[DrugRecommendation]:
    """
    Rank all compounds targeting a specific protein.

    Strategies: "efficacy" (IC50), "safety" (toxicity first), "balanced" (weight)
    """
    if target_id not in G:
        return []

    recommendations = []
    for compound_id in G.predecessors(target_id):
        node = G.nodes[compound_id]
        if node.get("node_type") != "compound":
            continue

        edge = G.edges[compound_id, target_id]
        ic50 = edge["ic50_nm"]
        source = edge.get("source", "experimental")
        phase = node.get("phase", 0)
        tox = 0.0
        if toxicity_data:
            tox = toxicity_data.get(compound_id, {}).get("overall", 0)
        else:
            tox = node.get("toxicity", 0)

        confidence, reasons = assign_confidence(ic50, source, phase, tox)

        recommendations.append(DrugRecommendation(
            rank=0,
            compound_id=compound_id,
            compound_name=node.get("name", compound_id),
            target_id=target_id,
            target_name=G.nodes[target_id].get("name", target_id),
            ic50_nm=ic50,
            weight=edge["weight"],
            toxicity=tox,
            source=source,
            confidence=confidence,
            confidence_reasons=reasons,
        ))

    # Compute QED and SA for multi-criteria sorting
    _qed_cache = {}
    _sa_cache = {}
    if strategy in ("safety", "balanced"):
        for rec in recommendations:
            smiles = G.nodes[rec.compound_id].get("smiles", "")
            if smiles:
                _qed_cache[rec.compound_id] = estimate_qed(smiles)
                _sa_cache[rec.compound_id] = estimate_synthetic_accessibility(smiles)

    sort_keys = {
        # Pure potency — IC50 only
        "efficacy": lambda r: r.ic50_nm,
        # Safety-first: penalize toxicity and poor drug-likeness, then IC50 as tiebreaker
        "safety": lambda r: (
            r.toxicity - _qed_cache.get(r.compound_id, 0.5) * 0.3,
            _sa_cache.get(r.compound_id, 5) / 10,
            r.ic50_nm,
        ),
        # Balanced: graph weight + drug-likeness bonus (lower QED = higher penalty)
        "balanced": lambda r: (
            r.weight - _qed_cache.get(r.compound_id, 0.5) * 0.5
            + _sa_cache.get(r.compound_id, 5) * 0.05
        ),
    }
    recommendations.sort(key=sort_keys.get(strategy, sort_keys["balanced"]))

    for i, rec in enumerate(recommendations):
        rec.rank = i + 1

    return recommendations


def compute_pareto_front(recommendations: list[DrugRecommendation]) -> list[DrugRecommendation]:
    """
    Pareto front over (IC50, toxicity). Both minimized.
    A solution is Pareto-optimal if no other is strictly better on both.
    """
    pareto = []
    for candidate in recommendations:
        dominated = False
        for other in recommendations:
            if other is candidate:
                continue
            if (other.ic50_nm <= candidate.ic50_nm and other.toxicity <= candidate.toxicity
                    and (other.ic50_nm < candidate.ic50_nm or other.toxicity < candidate.toxicity)):
                dominated = True
                break
        if not dominated:
            pareto.append(candidate)

    pareto.sort(key=lambda r: r.weight)
    for i, rec in enumerate(pareto):
        rec.rank = i + 1
    return pareto


def find_knee_point(pareto: list[DrugRecommendation]) -> DrugRecommendation | None:
    """Find the knee point — best balance between efficacy and toxicity."""
    if not pareto:
        return None
    if len(pareto) == 1:
        return pareto[0]

    ic50s = [r.ic50_nm for r in pareto]
    toxs = [r.toxicity for r in pareto]
    ic50_range = max(ic50s) - min(ic50s) or 1.0
    tox_range = max(toxs) - min(toxs) or 1.0

    best = None
    best_dist = float('inf')
    for r in pareto:
        norm_ic50 = (r.ic50_nm - min(ic50s)) / ic50_range
        norm_tox = (r.toxicity - min(toxs)) / tox_range
        dist = math.sqrt(norm_ic50 ** 2 + norm_tox ** 2)
        if dist < best_dist:
            best_dist = dist
            best = r
    return best


def compute_selectivity(
    G: nx.DiGraph,
    compound_id: str,
    primary_target: str,
    off_targets: list[str] | None = None,
) -> dict:
    """
    Compute selectivity profile for a compound.

    Selectivity ratio = IC50(off-target) / IC50(primary-target).
    Higher ratio = more selective for the primary target.
    """
    if compound_id not in G or primary_target not in G:
        return {"error": "Compound or target not found"}

    if not G.has_edge(compound_id, primary_target):
        return {"error": f"No activity data for {compound_id} vs {primary_target}"}

    primary_ic50 = G.edges[compound_id, primary_target]["ic50_nm"]

    # Find all targets this compound hits
    all_targets = []
    for target_id in G.successors(compound_id):
        node = G.nodes[target_id]
        if node.get("node_type") != "target":
            continue
        if target_id == primary_target:
            continue
        if off_targets and target_id not in off_targets:
            continue
        edge = G.edges[compound_id, target_id]
        ratio = edge["ic50_nm"] / primary_ic50 if primary_ic50 > 0 else float("inf")
        all_targets.append({
            "target_id": target_id,
            "target_name": node.get("name", target_id),
            "ic50_nm": edge["ic50_nm"],
            "selectivity_ratio": round(ratio, 1),
        })

    all_targets.sort(key=lambda x: x["selectivity_ratio"])

    return {
        "compound_id": compound_id,
        "compound_name": G.nodes[compound_id].get("name", compound_id),
        "primary_target": primary_target,
        "primary_ic50_nm": primary_ic50,
        "off_target_count": len(all_targets),
        "off_targets": all_targets,
        "is_selective": all(t["selectivity_ratio"] > 10 for t in all_targets) if all_targets else True,
    }


def compare_compounds_across_targets(
    G: nx.DiGraph,
    compound_ids: list[str],
    target_ids: list[str],
) -> list[dict]:
    """
    Build an IC50 matrix for comparing compounds across targets.
    Useful for selectivity profiling and polypharmacology analysis.
    """
    rows = []
    for cid in compound_ids:
        if cid not in G:
            continue
        row = {
            "compound_id": cid,
            "compound_name": G.nodes[cid].get("name", cid),
            "targets": {},
        }
        for tid in target_ids:
            if G.has_edge(cid, tid):
                edge = G.edges[cid, tid]
                row["targets"][G.nodes[tid].get("name", tid)] = edge["ic50_nm"]
            else:
                row["targets"][G.nodes[tid].get("name", tid)] = None
        rows.append(row)
    return rows


def format_recommendations(recommendations: list[DrugRecommendation], title: str = "Results") -> str:
    """Format recommendations as a readable report."""
    lines = [f"\n{'=' * 70}", f" {title}", f"{'=' * 70}"]
    for rec in recommendations:
        conf_badge = {"HIGH": "[+++]", "MEDIUM": "[++ ]", "LOW": "[+  ]"}[rec.confidence]
        if rec.ic50_nm < 0.1:
            ic50_str = f"{rec.ic50_nm:.3f}"
        elif rec.ic50_nm < 10:
            ic50_str = f"{rec.ic50_nm:.2f}"
        else:
            ic50_str = f"{rec.ic50_nm:.1f}"
        lines.append(
            f"  #{rec.rank} {conf_badge} {rec.compound_name:15s} -> {rec.target_name:6s} | "
            f"IC50={ic50_str:>8s}nM | Tox={rec.toxicity:.2f} | "
            f"Weight={rec.weight:.4f} | {rec.source:12s} | Conf={rec.confidence}"
        )
    lines.append(f"{'=' * 70}")
    return "\n".join(lines)


def format_confidence_detail(rec: DrugRecommendation) -> str:
    """Format detailed confidence explanation for a single recommendation."""
    lines = [f"\n  Confidence breakdown for {rec.compound_name} ({rec.confidence}):"]
    for reason in rec.confidence_reasons:
        lines.append(f"    - {reason}")
    return "\n".join(lines)
