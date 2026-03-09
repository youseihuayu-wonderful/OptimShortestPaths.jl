"""
Multi-hop biological network with probability-based edge weights.

Architecture:
  The network has 4 node layers and 4 edge types:

  Nodes:  compound → target → pathway → disease
  PPIs:   target ↔ target

Edge weight semantics:
  Every edge weight is a PROBABILITY in (0, 1]:
    - P(compound effectively binds target)    [from IC50]
    - P(target functionally drives pathway)   [from annotation confidence]
    - P(pathway causally linked to disease)   [from genetic/clinical evidence]
    - P(PPI interaction is real/functional)    [from STRING-like confidence]

  In log-space, w = -log(p), so:
    - Probabilities multiply along a path:  P(path) = product(p_i)
    - Log-costs add along a path:           w(path) = sum(w_i)
    - Shortest path = maximum probability path (most likely mechanism)

Multi-dimensional probability vectors:
  Each edge can carry a VECTOR of probabilities across independent dimensions:
    dim 0: efficacy   — P(this link contributes to therapeutic effect)
    dim 1: safety     — P(this link does NOT cause adverse effects)
    dim 2: evidence   — P(this link is real, based on data quality/source)

  Multi-objective shortest path finds Pareto-optimal paths across all dimensions.
"""

import math
import networkx as nx
from dataclasses import dataclass, field

from chempath.data.biological_network import get_biological_network


# ---------------------------------------------------------------------------
# Probability dimensions
# ---------------------------------------------------------------------------

PROB_DIMENSIONS = ["efficacy", "safety", "evidence"]


@dataclass
class EdgeProbability:
    """Probability vector for an edge. Each dimension is in (0, 1]."""
    efficacy: float = 1.0
    safety: float = 1.0
    evidence: float = 1.0

    def to_logcost(self) -> dict[str, float]:
        """Convert probabilities to additive log-costs: w = -log(p)."""
        return {
            "w_efficacy": -math.log(max(1e-10, self.efficacy)),
            "w_safety": -math.log(max(1e-10, self.safety)),
            "w_evidence": -math.log(max(1e-10, self.evidence)),
        }

    def to_dict(self) -> dict[str, float]:
        return {
            "p_efficacy": self.efficacy,
            "p_safety": self.safety,
            "p_evidence": self.evidence,
        }


# ---------------------------------------------------------------------------
# Probability transforms (domain data → uniform probability space)
# ---------------------------------------------------------------------------

def ic50_to_efficacy_prob(ic50_nm: float) -> float:
    """
    Convert IC50 (nM) to P(effective binding) via sigmoid.

    P = 1 / (1 + IC50 / IC50_ref)
    where IC50_ref = 100 nM (a typical "good" potency threshold).

    IC50=1nM    → P=0.99  (very likely to bind effectively)
    IC50=100nM  → P=0.50  (coin flip)
    IC50=10000nM→ P=0.01  (unlikely to be effective)
    """
    if ic50_nm <= 0:
        return 1e-10
    ic50_ref = 100.0  # nM
    return 1.0 / (1.0 + ic50_nm / ic50_ref)


def toxicity_to_safety_prob(tox_score: float) -> float:
    """
    Convert toxicity risk score [0,1] to P(safe).
    P(safe) = 1 - tox_score, clamped to [0.01, 1.0].
    """
    return max(0.01, min(1.0, 1.0 - tox_score))


def phase_to_evidence_prob(phase: int, source: str = "experimental") -> float:
    """
    Convert clinical phase + data source to P(evidence is reliable).

    Phase 4 (approved) + experimental  → 0.95
    Phase 3 + experimental              → 0.80
    Phase 1-2 + experimental            → 0.60
    Preclinical + experimental           → 0.40
    Predicted                            → 0.20
    """
    if source != "experimental":
        return 0.20

    phase_probs = {0: 0.40, 1: 0.55, 2: 0.65, 3: 0.80, 4: 0.95}
    return phase_probs.get(min(phase, 4), 0.40)


def annotation_to_prob(weight_01: float) -> float:
    """
    Convert a curated annotation weight [0, 1] to probability.
    Input weight is COST (lower = stronger); we invert: P = 1 - weight.
    Clamped to [0.05, 0.99].
    """
    return max(0.05, min(0.99, 1.0 - weight_01))


# ---------------------------------------------------------------------------
# Multi-hop graph construction
# ---------------------------------------------------------------------------

def build_multihop_graph(
    chembl_data: dict,
    toxicity_penalty: float = 1.0,
    verbose: bool = True,
) -> nx.DiGraph:
    """
    Build a 4-layer directed graph:
      compound → target → pathway → disease
    with PPI edges between targets.

    Every edge carries a 3-dimensional probability vector:
      (p_efficacy, p_safety, p_evidence)
    and corresponding log-costs:
      (w_efficacy, w_safety, w_evidence)
    """
    G = nx.DiGraph()
    bio = get_biological_network()

    # ------- Layer 1: Compounds -------
    compounds_by_id = {}
    for compound in chembl_data.get("compounds", []):
        cid = compound["chembl_id"]
        compounds_by_id[cid] = compound
        G.add_node(cid,
                   name=compound.get("name", cid),
                   node_type="compound",
                   layer=0,
                   smiles=compound.get("smiles", ""),
                   phase=compound.get("phase", 0))

    # ------- Layer 2: Targets -------
    target_ids = set()
    for target in chembl_data.get("targets", []):
        tid = target["chembl_id"]
        target_ids.add(tid)
        G.add_node(tid,
                   name=target.get("name", tid),
                   node_type="target",
                   layer=1,
                   organism=target.get("organism", ""))

    # ------- Layer 3: Pathways -------
    for pathway in bio["pathways"]:
        G.add_node(pathway["id"],
                   name=pathway["name"],
                   node_type="pathway",
                   layer=2,
                   description=pathway.get("description", ""),
                   kegg_id=pathway.get("kegg_id", ""))

    # ------- Layer 4: Diseases -------
    for disease in bio["diseases"]:
        G.add_node(disease["id"],
                   name=disease["name"],
                   node_type="disease",
                   layer=3,
                   full_name=disease.get("full_name", ""),
                   icd10=disease.get("icd10", ""))

    # ------- Edges: Compound → Target (from bioactivity data) -------
    toxicity_data = chembl_data.get("toxicity", {})
    edge_count = {"compound_target": 0, "target_pathway": 0,
                  "pathway_disease": 0, "ppi": 0}

    for activity in chembl_data.get("bioactivities", []):
        cid = activity["compound"]
        tid = activity["target"]
        if cid not in G or tid not in G:
            continue

        ic50 = activity["value"]
        source = activity.get("source", "experimental")
        phase = compounds_by_id.get(cid, {}).get("phase", 0)
        tox = toxicity_data.get(cid, {}).get("overall", 0)

        prob = EdgeProbability(
            efficacy=ic50_to_efficacy_prob(ic50),
            safety=toxicity_to_safety_prob(tox),
            evidence=phase_to_evidence_prob(phase, source),
        )
        logcost = prob.to_logcost()

        G.add_edge(cid, tid,
                   edge_type="compound_target",
                   weight=logcost["w_efficacy"],
                   ic50_nm=ic50,
                   **prob.to_dict(),
                   **logcost,
                   source=source)
        edge_count["compound_target"] += 1

    # ------- Edges: Target → Pathway -------
    for edge in bio["target_pathway_edges"]:
        src, tgt = edge["source"], edge["target"]
        if src not in G or tgt not in G:
            continue

        p_eff = annotation_to_prob(edge["weight"])
        prob = EdgeProbability(
            efficacy=p_eff,
            safety=0.95,  # pathway membership doesn't affect safety much
            evidence=0.85,  # curated from KEGG, high evidence
        )
        logcost = prob.to_logcost()

        G.add_edge(src, tgt,
                   edge_type="target_pathway",
                   weight=logcost["w_efficacy"],
                   role=edge.get("role", ""),
                   biological_evidence=edge.get("evidence", ""),
                   **prob.to_dict(),
                   **logcost)
        edge_count["target_pathway"] += 1

    # ------- Edges: Pathway → Disease -------
    for edge in bio["pathway_disease_edges"]:
        src, tgt = edge["source"], edge["target"]
        if src not in G or tgt not in G:
            continue

        p_eff = annotation_to_prob(edge["weight"])
        prob = EdgeProbability(
            efficacy=p_eff,
            safety=0.90,  # pathway-disease link doesn't directly affect safety
            evidence=0.80,  # curated from KEGG/DisGeNET
        )
        logcost = prob.to_logcost()

        G.add_edge(src, tgt,
                   edge_type="pathway_disease",
                   weight=logcost["w_efficacy"],
                   biological_evidence=edge.get("evidence", ""),
                   **prob.to_dict(),
                   **logcost)
        edge_count["pathway_disease"] += 1

    # ------- Edges: Target ↔ Target (PPI) -------
    for edge in bio["ppi_edges"]:
        src, tgt = edge["source"], edge["target"]
        if src not in G or tgt not in G:
            continue

        p_eff = annotation_to_prob(edge["weight"])
        ppi_type = edge.get("type", "interaction")

        # Safety depends on PPI type: cross-resistance is a safety concern
        p_safe = 0.70 if ppi_type == "cross_resistance" else 0.90
        # Negative feedback loops reduce evidence certainty
        p_evid = 0.60 if ppi_type == "negative_feedback" else 0.75

        prob = EdgeProbability(efficacy=p_eff, safety=p_safe, evidence=p_evid)
        logcost = prob.to_logcost()

        G.add_edge(src, tgt,
                   edge_type="ppi",
                   weight=logcost["w_efficacy"],
                   ppi_type=ppi_type,
                   biological_evidence=edge.get("evidence", ""),
                   **prob.to_dict(),
                   **logcost)
        edge_count["ppi"] += 1

        # Add reverse edge for bidirectional PPIs
        if edge.get("bidirectional", False):
            G.add_edge(tgt, src,
                       edge_type="ppi",
                       weight=logcost["w_efficacy"],
                       ppi_type=ppi_type,
                       biological_evidence=edge.get("evidence", ""),
                       **prob.to_dict(),
                       **logcost)
            edge_count["ppi"] += 1

    if verbose:
        n = G.number_of_nodes()
        e = G.number_of_edges()
        print(f"  Multi-hop graph: {n} nodes, {e} edges")
        for k, v in edge_count.items():
            print(f"    {k}: {v} edges")

    return G


def get_multihop_summary(G: nx.DiGraph) -> dict:
    """Summary statistics for the multi-hop graph."""
    layers = {"compound": 0, "target": 0, "pathway": 0, "disease": 0}
    for _, d in G.nodes(data=True):
        nt = d.get("node_type", "unknown")
        if nt in layers:
            layers[nt] += 1

    edge_types = {}
    for _, _, d in G.edges(data=True):
        et = d.get("edge_type", "unknown")
        edge_types[et] = edge_types.get(et, 0) + 1

    return {
        "total_nodes": G.number_of_nodes(),
        "total_edges": G.number_of_edges(),
        "nodes_by_type": layers,
        "edges_by_type": edge_types,
    }
