"""
Graph construction from drug-target bioactivity data.

Transformation:
  - Compounds and targets → nodes
  - Bioactivities → edges
  - Edge weight = inverted pIC50 + toxicity penalty (lower = better)
"""

import math
import networkx as nx
from chempath.chemistry.smiles import filter_valid_compounds
from chempath.chemistry.properties import compute_risk_scores


def ic50_to_weight(ic50_nm: float) -> float:
    """
    Convert IC50 (nM) to edge weight for shortest-path optimization.

    Uses -log10 transform: weight = -log10(IC50_nM * 1e-9) mapped so that
    lower IC50 (higher potency) → lower weight.

    Scale: IC50=0.01nM → ~1.0, IC50=1nM → ~1.7, IC50=100nM → ~3.0, IC50=10000nM → ~5.0
    """
    if ic50_nm <= 0:
        return float('inf')
    pic50 = -math.log10(ic50_nm * 1e-9)  # pIC50: higher = more potent
    # Map pIC50 to weight: weight = 14 - pIC50, clamped to [0.01, 10.0]
    # This ensures sub-picomolar compounds still get distinct (small) weights
    weight = max(0.01, min(10.0, 14.0 - pic50))
    return round(weight, 4)


def build_drug_target_graph(
    data: dict,
    toxicity_penalty: float = 1.0,
    verbose: bool = True,
) -> nx.DiGraph:
    """Build a directed graph from bioactivity data with SMILES validation."""
    G = nx.DiGraph()

    valid_compounds, rejected = filter_valid_compounds(data["compounds"])
    if rejected and verbose:
        print(f"  [SMILES Validator] Rejected {len(rejected)} compounds with invalid SMILES:")
        for r in rejected:
            print(f"    - {r['name']}: {r['_smiles_validation']}")

    # Use explicit toxicity data if provided, otherwise compute risk scores from SMILES
    toxicity_data = data.get("toxicity", {})
    if not toxicity_data:
        risk_scores = compute_risk_scores(valid_compounds)
        toxicity_data = {cid: {"overall": score} for cid, score in risk_scores.items()}

    for compound in valid_compounds:
        cid = compound["chembl_id"]
        G.add_node(
            cid,
            name=compound["name"],
            node_type="compound",
            smiles=compound["smiles"],
            phase=compound.get("phase", 0),
            toxicity=toxicity_data.get(cid, {}).get("overall", 0),
        )

    for target in data["targets"]:
        G.add_node(
            target["chembl_id"],
            name=target["name"],
            node_type="target",
            organism=target.get("organism", ""),
        )

    for activity in data["bioactivities"]:
        compound_id = activity["compound"]
        target_id = activity["target"]

        if compound_id not in G:
            continue

        ic50 = activity["value"]
        base_weight = ic50_to_weight(ic50)
        tox = toxicity_data.get(compound_id, {}).get("overall", 0)
        weight = base_weight + tox * toxicity_penalty

        G.add_edge(
            compound_id, target_id,
            weight=round(weight, 4),
            base_weight=base_weight,
            ic50_nm=ic50,
            toxicity_penalty=round(tox * toxicity_penalty, 4),
            source=activity.get("source", "experimental"),
        )

    return G


def add_predicted_edges(
    G: nx.DiGraph,
    predictions: list[dict],
    uncertainty_penalty: float = 0.2,
) -> nx.DiGraph:
    """Add Chemprop-predicted edges with an uncertainty penalty."""
    for pred in predictions:
        compound_id = pred["compound"]
        target_id = pred["target"]

        if compound_id not in G or target_id not in G:
            continue
        if G.has_edge(compound_id, target_id):
            continue

        ic50 = pred["predicted_ic50"]
        base_weight = ic50_to_weight(ic50)
        weight = base_weight + uncertainty_penalty

        G.add_edge(
            compound_id, target_id,
            weight=round(weight, 4),
            base_weight=base_weight,
            ic50_nm=ic50,
            uncertainty_penalty=uncertainty_penalty,
            source="predicted",
        )

    return G


def get_graph_summary(G: nx.DiGraph) -> dict:
    """Return a summary of the graph structure."""
    compounds = [n for n, d in G.nodes(data=True) if d.get("node_type") == "compound"]
    targets = [n for n, d in G.nodes(data=True) if d.get("node_type") == "target"]
    exp_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get("source") == "experimental"]
    pred_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get("source") == "predicted"]

    return {
        "total_nodes": G.number_of_nodes(),
        "compounds": len(compounds),
        "targets": len(targets),
        "total_edges": G.number_of_edges(),
        "experimental_edges": len(exp_edges),
        "predicted_edges": len(pred_edges),
    }
