"""
Multi-hop pathfinding algorithms for drug discovery.

Four core use cases where shortest-path is genuinely the right algorithm:

1. Drug Repurposing:   compound → target → pathway → disease
   "Can an existing drug treat a new disease via an indirect mechanism?"

2. Mechanism of Action: compound → ??? → ??? → disease
   "What is the most likely biological mechanism for this drug-disease link?"

3. Target Identification: disease → ??? → target
   "Which protein target is most reachable from this disease phenotype?"

4. Combination Therapy: find minimum-cost set of paths covering {T1, T2, T3}
   "What drug combination optimally covers multiple targets?"

All path costs are in -log(probability) space, so:
  - Shorter path = higher probability = more likely mechanism
  - Path cost = -log(P(drug reaches disease through this mechanism))
"""

import math
import heapq
import networkx as nx
from dataclasses import dataclass, field


@dataclass
class PathResult:
    """A single path through the biological network."""
    path: list[str]            # node IDs in order
    path_names: list[str]      # human-readable names
    path_types: list[str]      # node types (compound, target, pathway, disease)
    edge_types: list[str]      # edge types along the path
    total_cost: float          # sum of -log(p) = -log(product of probabilities)
    total_probability: float   # product of probabilities along path
    costs_by_dim: dict[str, float]  # per-dimension costs
    probs_by_dim: dict[str, float]  # per-dimension probabilities
    edges_detail: list[dict]   # detailed info per edge

    @property
    def n_hops(self) -> int:
        return len(self.path) - 1


@dataclass
class RepurposingResult:
    """Result of a drug repurposing query."""
    compound_id: str
    compound_name: str
    disease_id: str
    disease_name: str
    paths: list[PathResult]
    best_probability: float
    pareto_paths: list[PathResult]  # Pareto-optimal across dimensions


@dataclass
class MoAResult:
    """Result of a mechanism-of-action query."""
    compound_id: str
    compound_name: str
    disease_id: str
    disease_name: str
    paths: list[PathResult]      # all shortest paths (may include ties)
    primary_mechanism: PathResult | None  # single best path


@dataclass
class TargetResult:
    """Result of a target identification query."""
    disease_id: str
    disease_name: str
    targets: list[dict]   # [{target_id, target_name, best_path, probability}, ...]


@dataclass
class CombinationResult:
    """Result of a combination therapy query."""
    target_ids: list[str]
    target_names: list[str]
    compounds: list[dict]  # [{compound_id, compound_name, targets_covered, paths, total_cost}]
    best_combination: list[dict]  # minimum set of compounds covering all targets


# ---------------------------------------------------------------------------
# Core pathfinding
# ---------------------------------------------------------------------------

def _reconstruct_path(G: nx.DiGraph, predecessors: dict, source: str, target: str) -> list[str]:
    """Reconstruct path from predecessor map."""
    path = []
    node = target
    while node is not None:
        path.append(node)
        node = predecessors.get(node)
    path.reverse()
    if path and path[0] == source:
        return path
    return []


def _build_path_result(G: nx.DiGraph, path: list[str], weight_key: str = "w_efficacy") -> PathResult:
    """Build a PathResult from a list of node IDs."""
    path_names = [G.nodes[n].get("name", n) for n in path]
    path_types = [G.nodes[n].get("node_type", "unknown") for n in path]

    edge_types = []
    edges_detail = []
    costs = {"efficacy": 0.0, "safety": 0.0, "evidence": 0.0}

    for i in range(len(path) - 1):
        u, v = path[i], path[i + 1]
        edge_data = G.edges[u, v]
        edge_types.append(edge_data.get("edge_type", "unknown"))

        w_eff = edge_data.get("w_efficacy", 0)
        w_saf = edge_data.get("w_safety", 0)
        w_evid = edge_data.get("w_evidence", 0)
        costs["efficacy"] += w_eff
        costs["safety"] += w_saf
        costs["evidence"] += w_evid

        edges_detail.append({
            "from": G.nodes[u].get("name", u),
            "to": G.nodes[v].get("name", v),
            "from_id": u,
            "to_id": v,
            "edge_type": edge_data.get("edge_type", ""),
            "p_efficacy": edge_data.get("p_efficacy", 1.0),
            "p_safety": edge_data.get("p_safety", 1.0),
            "p_evidence": edge_data.get("p_evidence", 1.0),
            "ic50_nm": edge_data.get("ic50_nm"),
            "biological_evidence": edge_data.get("biological_evidence", ""),
        })

    probs = {dim: math.exp(-cost) for dim, cost in costs.items()}
    total_cost = costs.get(weight_key.replace("w_", ""), costs["efficacy"])
    total_prob = math.exp(-total_cost)

    return PathResult(
        path=path,
        path_names=path_names,
        path_types=path_types,
        edge_types=edge_types,
        total_cost=total_cost,
        total_probability=total_prob,
        costs_by_dim=costs,
        probs_by_dim=probs,
        edges_detail=edges_detail,
    )


def find_shortest_paths(
    G: nx.DiGraph,
    source: str,
    target: str,
    weight_key: str = "w_efficacy",
    k: int = 5,
) -> list[PathResult]:
    """
    Find top-k shortest paths from source to target using Yen's algorithm.

    weight_key: which probability dimension to optimize
      "w_efficacy"  — maximize P(therapeutic effect along path)
      "w_safety"    — maximize P(no adverse effects along path)
      "w_evidence"  — maximize P(all links are real)
    """
    if source not in G or target not in G:
        return []

    results = []

    # Find first shortest path with Dijkstra
    try:
        first_path = nx.shortest_path(G, source, target, weight=weight_key)
    except (nx.NetworkXNoPath, nx.NodeNotFound):
        return []

    results.append(_build_path_result(G, first_path, weight_key))

    # Yen's k-shortest paths
    if k > 1:
        A = [first_path]  # accepted paths
        B = []            # candidate paths (min-heap)

        for ki in range(1, k):
            for i in range(len(A[-1]) - 1):
                spur_node = A[-1][i]
                root_path = A[-1][:i + 1]
                root_cost = sum(
                    G.edges[root_path[j], root_path[j + 1]].get(weight_key, float("inf"))
                    for j in range(len(root_path) - 1)
                )

                # Remove edges that share the same root path
                removed_edges = []
                for path in A:
                    if path[:i + 1] == root_path and i + 1 < len(path):
                        u, v = path[i], path[i + 1]
                        if G.has_edge(u, v):
                            removed_edges.append((u, v, G.edges[u, v].copy()))
                            G.remove_edge(u, v)

                # Remove root path nodes (except spur node)
                removed_nodes = []
                for node in root_path[:-1]:
                    if node in G and node != spur_node:
                        removed_nodes.append((node, dict(G.nodes[node]),
                                              list(G.in_edges(node, data=True)),
                                              list(G.out_edges(node, data=True))))
                        G.remove_node(node)

                try:
                    spur_path = nx.shortest_path(G, spur_node, target, weight=weight_key)
                    total_path = root_path[:-1] + spur_path
                    total_cost = root_cost + sum(
                        G.edges[spur_path[j], spur_path[j + 1]].get(weight_key, float("inf"))
                        for j in range(len(spur_path) - 1)
                    )
                    if total_path not in A:
                        heapq.heappush(B, (total_cost, total_path))
                except (nx.NetworkXNoPath, nx.NodeNotFound):
                    pass

                # Restore removed edges and nodes
                for node, attrs, in_edges, out_edges in reversed(removed_nodes):
                    G.add_node(node, **attrs)
                    for u, v, d in in_edges:
                        if u in G:
                            G.add_edge(u, v, **d)
                    for u, v, d in out_edges:
                        if v in G:
                            G.add_edge(u, v, **d)
                for u, v, d in removed_edges:
                    G.add_edge(u, v, **d)

            if not B:
                break

            _, next_path = heapq.heappop(B)
            A.append(next_path)
            results.append(_build_path_result(G, next_path, weight_key))

    return results


# ---------------------------------------------------------------------------
# Use Case 1: Drug Repurposing
# ---------------------------------------------------------------------------

def find_repurposing_candidates(
    G: nx.DiGraph,
    compound_id: str,
    k_paths: int = 5,
    weight_key: str = "w_efficacy",
) -> list[RepurposingResult]:
    """
    Find all diseases reachable from a compound through the biological network.
    Returns ranked list of (compound, disease) pairs with best paths.

    This is the core drug repurposing query: "What diseases could this drug treat?"
    """
    if compound_id not in G:
        return []

    compound_name = G.nodes[compound_id].get("name", compound_id)

    # Find all disease nodes
    disease_nodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "disease"]

    results = []
    for disease_id in disease_nodes:
        paths = find_shortest_paths(G, compound_id, disease_id,
                                    weight_key=weight_key, k=k_paths)
        if not paths:
            continue

        # Compute Pareto-optimal paths across dimensions
        pareto = _compute_pareto_paths(paths)

        results.append(RepurposingResult(
            compound_id=compound_id,
            compound_name=compound_name,
            disease_id=disease_id,
            disease_name=G.nodes[disease_id].get("name", disease_id),
            paths=paths,
            best_probability=paths[0].total_probability,
            pareto_paths=pareto,
        ))

    # Sort by best probability (highest first)
    results.sort(key=lambda r: -r.best_probability)
    return results


# ---------------------------------------------------------------------------
# Use Case 2: Mechanism of Action
# ---------------------------------------------------------------------------

def find_mechanism_of_action(
    G: nx.DiGraph,
    compound_id: str,
    disease_id: str,
    k_paths: int = 10,
) -> MoAResult | None:
    """
    Find the most likely biological mechanism(s) connecting a drug to a disease.

    Returns all shortest paths (candidate mechanisms), with the primary mechanism
    being the path with the highest combined probability.
    """
    if compound_id not in G or disease_id not in G:
        return None

    # Find paths optimizing different dimensions
    efficacy_paths = find_shortest_paths(G, compound_id, disease_id,
                                         weight_key="w_efficacy", k=k_paths)
    safety_paths = find_shortest_paths(G, compound_id, disease_id,
                                        weight_key="w_safety", k=k_paths)
    evidence_paths = find_shortest_paths(G, compound_id, disease_id,
                                          weight_key="w_evidence", k=k_paths)

    # Merge and deduplicate paths
    all_paths = {}
    for path in efficacy_paths + safety_paths + evidence_paths:
        key = tuple(path.path)
        if key not in all_paths:
            all_paths[key] = path

    paths = sorted(all_paths.values(), key=lambda p: p.total_cost)

    if not paths:
        return None

    # Primary mechanism: best combined probability (geometric mean across dims)
    best = min(paths, key=lambda p: sum(p.costs_by_dim.values()))

    return MoAResult(
        compound_id=compound_id,
        compound_name=G.nodes[compound_id].get("name", compound_id),
        disease_id=disease_id,
        disease_name=G.nodes[disease_id].get("name", disease_id),
        paths=paths,
        primary_mechanism=best,
    )


# ---------------------------------------------------------------------------
# Use Case 3: Target Identification
# ---------------------------------------------------------------------------

def identify_targets_for_disease(
    G: nx.DiGraph,
    disease_id: str,
    weight_key: str = "w_efficacy",
) -> TargetResult | None:
    """
    Find the most reachable protein targets from a disease (reverse direction).

    Uses reverse graph: disease → pathway → target.
    Identifies which targets are most strongly connected to the disease.
    """
    if disease_id not in G:
        return None

    # Build reverse graph for backward search
    G_rev = G.reverse()

    target_nodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "target"]

    targets = []
    for tid in target_nodes:
        try:
            path = nx.shortest_path(G_rev, disease_id, tid, weight=weight_key)
            # Convert reverse path to forward path
            fwd_path = list(reversed(path))
            result = _build_path_result(G, fwd_path, weight_key)
            targets.append({
                "target_id": tid,
                "target_name": G.nodes[tid].get("name", tid),
                "best_path": result,
                "probability": result.total_probability,
                "n_hops": result.n_hops,
                "mechanism": " → ".join(result.path_names),
            })
        except (nx.NetworkXNoPath, nx.NodeNotFound):
            continue

    targets.sort(key=lambda t: -t["probability"])

    return TargetResult(
        disease_id=disease_id,
        disease_name=G.nodes[disease_id].get("name", disease_id),
        targets=targets,
    )


# ---------------------------------------------------------------------------
# Use Case 4: Combination Therapy
# ---------------------------------------------------------------------------

def find_combination_therapy(
    G: nx.DiGraph,
    target_ids: list[str],
    disease_id: str | None = None,
    max_compounds: int = 3,
    weight_key: str = "w_efficacy",
) -> CombinationResult | None:
    """
    Find the minimum set of compounds that covers all specified targets.

    This is a set cover / Steiner tree problem:
    - For each target, find compounds that reach it
    - Find the minimum set of compounds that covers all targets
    - Optionally require paths to also reach a specific disease
    """
    valid_targets = [t for t in target_ids if t in G]
    if not valid_targets:
        return None

    target_names = [G.nodes[t].get("name", t) for t in valid_targets]

    # Find compounds that reach each target
    compound_nodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "compound"]

    compound_coverage = {}
    for cid in compound_nodes:
        covered = []
        total_cost = 0.0
        paths = []

        for tid in valid_targets:
            if G.has_edge(cid, tid):
                edge = G.edges[cid, tid]
                cost = edge.get(weight_key, float("inf"))
                covered.append(tid)
                total_cost += cost

                # If disease specified, find full path
                if disease_id and disease_id in G:
                    full_paths = find_shortest_paths(G, cid, disease_id,
                                                     weight_key=weight_key, k=1)
                    if full_paths:
                        paths.append(full_paths[0])
                else:
                    # Build single-hop path result
                    path_result = _build_path_result(G, [cid, tid], weight_key)
                    paths.append(path_result)

        if covered:
            compound_coverage[cid] = {
                "compound_id": cid,
                "compound_name": G.nodes[cid].get("name", cid),
                "targets_covered": covered,
                "target_names_covered": [G.nodes[t].get("name", t) for t in covered],
                "n_covered": len(covered),
                "total_cost": total_cost,
                "avg_cost": total_cost / len(covered),
                "paths": paths,
            }

    # Greedy set cover: find minimum compounds covering all targets
    uncovered = set(valid_targets)
    combination = []
    used = set()

    while uncovered and len(combination) < max_compounds:
        # Pick compound covering the most uncovered targets, breaking ties by cost
        best = None
        best_score = (-1, float("inf"))

        for cid, info in compound_coverage.items():
            if cid in used:
                continue
            new_covered = len(set(info["targets_covered"]) & uncovered)
            if new_covered == 0:
                continue
            score = (new_covered, -info["avg_cost"])
            if score > best_score:
                best_score = score
                best = cid

        if best is None:
            break

        combination.append(compound_coverage[best])
        uncovered -= set(compound_coverage[best]["targets_covered"])
        used.add(best)

    # Sort all compounds by coverage then cost
    all_compounds = sorted(compound_coverage.values(),
                           key=lambda c: (-c["n_covered"], c["avg_cost"]))

    return CombinationResult(
        target_ids=valid_targets,
        target_names=target_names,
        compounds=all_compounds[:20],
        best_combination=combination,
    )


# ---------------------------------------------------------------------------
# Pareto path computation
# ---------------------------------------------------------------------------

def _compute_pareto_paths(paths: list[PathResult]) -> list[PathResult]:
    """Find Pareto-optimal paths across the 3 probability dimensions."""
    pareto = []
    for candidate in paths:
        dominated = False
        for other in paths:
            if other is candidate:
                continue
            # Check if other dominates candidate (lower cost in all dimensions)
            if (other.costs_by_dim["efficacy"] <= candidate.costs_by_dim["efficacy"] and
                other.costs_by_dim["safety"] <= candidate.costs_by_dim["safety"] and
                other.costs_by_dim["evidence"] <= candidate.costs_by_dim["evidence"] and
                (other.costs_by_dim["efficacy"] < candidate.costs_by_dim["efficacy"] or
                 other.costs_by_dim["safety"] < candidate.costs_by_dim["safety"] or
                 other.costs_by_dim["evidence"] < candidate.costs_by_dim["evidence"])):
                dominated = True
                break
        if not dominated:
            pareto.append(candidate)
    return pareto


# ---------------------------------------------------------------------------
# Convenience: scan all compound-disease pairs
# ---------------------------------------------------------------------------

def scan_all_repurposing(
    G: nx.DiGraph,
    weight_key: str = "w_efficacy",
    min_probability: float = 0.001,
) -> list[dict]:
    """
    Scan all (compound, disease) pairs for repurposing opportunities.
    Returns a ranked list of the most promising connections.
    """
    compound_nodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "compound"]
    disease_nodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "disease"]

    opportunities = []

    for disease_id in disease_nodes:
        disease_name = G.nodes[disease_id].get("name", disease_id)

        for cid in compound_nodes:
            paths = find_shortest_paths(G, cid, disease_id,
                                         weight_key=weight_key, k=1)
            if not paths:
                continue

            best = paths[0]
            if best.total_probability < min_probability:
                continue

            opportunities.append({
                "compound_id": cid,
                "compound_name": G.nodes[cid].get("name", cid),
                "disease_id": disease_id,
                "disease_name": disease_name,
                "probability": best.total_probability,
                "n_hops": best.n_hops,
                "mechanism": " → ".join(best.path_names),
                "path": best,
                "probs_by_dim": best.probs_by_dim,
            })

    opportunities.sort(key=lambda o: -o["probability"])
    return opportunities
