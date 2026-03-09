"""
Bridge between ChemPath (Python/NetworkX) and OptimShortestPaths.jl (Julia/DMY).

Exports the multi-hop biological network as an edge list that the Julia DMY
algorithm can ingest, then runs shortest-path queries using the O(m log^(2/3) n)
solver and returns results back to Python.

Usage:
    from chempath.graph.julia_bridge import JuliaBridge

    bridge = JuliaBridge(G)            # G is a NetworkX DiGraph from build_multihop_graph()
    result = bridge.run_sssp("CHEMBL939", weight_key="w_efficacy")
    print(result.distances["DIS_NSCLC"])   # distance to NSCLC
    print(result.paths["DIS_NSCLC"])       # path as list of node IDs
"""

import json
import subprocess
import tempfile
import math
from pathlib import Path
from dataclasses import dataclass, field

import networkx as nx


# Path to the Julia project root (parent of ChemPath/)
JULIA_PROJECT = Path(__file__).resolve().parent.parent.parent.parent
JULIA_BRIDGE_SCRIPT = Path(__file__).resolve().parent / "dmy_bridge.jl"


def find_julia() -> str | None:
    """Find the Julia binary. Returns path or None."""
    import shutil
    path = shutil.which("julia")
    if path:
        return path
    for candidate in [
        Path.home() / ".juliaup" / "bin" / "julia",
        Path("/opt/homebrew/bin/julia"),
        Path("/usr/local/bin/julia"),
    ]:
        if candidate.exists():
            return str(candidate)
    return None


JULIA_BINARY = find_julia()


@dataclass
class DMYResult:
    """Result from Julia DMY shortest-path computation."""
    source: str
    weight_key: str
    distances: dict[str, float]    # node_id → distance
    parents: dict[str, str | None]  # node_id → parent_node_id
    paths: dict[str, list[str]]    # node_id → [source, ..., node_id]
    node_names: dict[str, str]     # node_id → human-readable name
    n_vertices: int
    n_edges: int
    n_reachable: int
    algorithm: str = "DMY O(m log^(2/3) n)"


class JuliaBridge:
    """Bridge to call OptimShortestPaths.jl DMY solver on ChemPath graphs."""

    def __init__(self, G: nx.DiGraph):
        self.G = G
        self._node_to_idx: dict[str, int] = {}
        self._idx_to_node: dict[int, str] = {}
        self._build_index()

    def _build_index(self):
        """Map NetworkX node IDs (strings) to Julia vertex indices (1-based ints)."""
        for i, node in enumerate(self.G.nodes(), start=1):
            self._node_to_idx[node] = i
            self._idx_to_node[i] = node

    def export_graph(self, weight_key: str = "w_efficacy") -> dict:
        """Export NetworkX graph as JSON-serializable edge list for Julia."""
        edges = []
        for u, v, d in self.G.edges(data=True):
            w = d.get(weight_key, d.get("weight", 0.0))
            edges.append({
                "source": self._node_to_idx[u],
                "target": self._node_to_idx[v],
                "weight": max(0.0, w),  # DMY requires non-negative
            })

        # Node metadata
        nodes = {}
        for node_id, idx in self._node_to_idx.items():
            nd = self.G.nodes[node_id]
            nodes[str(idx)] = {
                "id": node_id,
                "name": nd.get("name", node_id),
                "type": nd.get("node_type", "unknown"),
            }

        return {
            "n_vertices": len(self._node_to_idx),
            "edges": edges,
            "nodes": nodes,
        }

    def run_sssp(
        self,
        source: str,
        weight_key: str = "w_efficacy",
        targets: list[str] | None = None,
        timeout: int = 60,
    ) -> DMYResult:
        """
        Run DMY single-source shortest paths from source to all vertices.

        Args:
            source: source node ID (e.g. "CHEMBL939")
            weight_key: which weight dimension to use
            targets: if provided, only return distances/paths for these nodes
            timeout: seconds before killing Julia process

        Returns:
            DMYResult with distances, parents, and reconstructed paths
        """
        if JULIA_BINARY is None:
            raise RuntimeError(
                "Julia not found. Install Julia (https://julialang.org/downloads/) "
                "and ensure it is on PATH to use the DMY bridge."
            )

        if source not in self._node_to_idx:
            raise ValueError(f"Source node '{source}' not in graph")

        source_idx = self._node_to_idx[source]
        graph_data = self.export_graph(weight_key)
        graph_data["source"] = source_idx

        # Write input to temp file
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".json", delete=False, prefix="chempath_dmy_"
        ) as f:
            json.dump(graph_data, f)
            input_path = f.name

        output_path = input_path.replace(".json", "_result.json")

        try:
            # Call Julia
            cmd = [
                JULIA_BINARY,
                f"--project={JULIA_PROJECT}",
                str(JULIA_BRIDGE_SCRIPT),
                input_path,
                output_path,
            ]

            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout,
                cwd=str(JULIA_PROJECT),
            )

            if proc.returncode != 0:
                raise RuntimeError(
                    f"Julia DMY solver failed (exit {proc.returncode}):\n"
                    f"stdout: {proc.stdout[:500]}\n"
                    f"stderr: {proc.stderr[:500]}"
                )

            # Read results
            with open(output_path) as f:
                result = json.load(f)

            return self._parse_result(result, source, weight_key, targets)

        finally:
            # Clean up temp files
            Path(input_path).unlink(missing_ok=True)
            Path(output_path).unlink(missing_ok=True)

    def _parse_result(
        self,
        result: dict,
        source: str,
        weight_key: str,
        targets: list[str] | None,
    ) -> DMYResult:
        """Parse Julia output back into Python objects."""
        distances_raw = result["distances"]  # list indexed by vertex (1-based)
        parents_raw = result["parents"]      # list indexed by vertex (1-based)

        distances = {}
        parents = {}
        paths = {}
        node_names = {}
        n_reachable = 0

        for idx_str, dist in enumerate(distances_raw, start=1):
            node_id = self._idx_to_node.get(idx_str)
            if node_id is None:
                continue

            if targets and node_id not in targets and node_id != source:
                continue

            node_names[node_id] = self.G.nodes[node_id].get("name", node_id)

            if dist < 1e300:  # reachable (Julia INF = typemax(Float64))
                distances[node_id] = dist
                n_reachable += 1

                # Reconstruct path from parents
                path = self._reconstruct_path(parents_raw, idx_str, self._node_to_idx[source])
                paths[node_id] = [self._idx_to_node[i] for i in path]

                parent_idx = parents_raw[idx_str - 1]  # 0-indexed in list
                parents[node_id] = self._idx_to_node.get(parent_idx) if parent_idx > 0 else None
            else:
                distances[node_id] = float("inf")
                parents[node_id] = None
                paths[node_id] = []

        return DMYResult(
            source=source,
            weight_key=weight_key,
            distances=distances,
            parents=parents,
            paths=paths,
            node_names=node_names,
            n_vertices=result.get("n_vertices", len(self._node_to_idx)),
            n_edges=result.get("n_edges", self.G.number_of_edges()),
            n_reachable=n_reachable,
        )

    def _reconstruct_path(
        self, parents: list[int], target_idx: int, source_idx: int
    ) -> list[int]:
        """Reconstruct path from parent array (1-based indices)."""
        path = []
        current = target_idx
        visited = set()
        while current > 0 and current not in visited:
            visited.add(current)
            path.append(current)
            if current == source_idx:
                break
            parent = parents[current - 1]  # parents list is 0-indexed
            current = parent
        path.reverse()
        if path and path[0] == source_idx:
            return path
        return [source_idx, target_idx]  # fallback: direct edge

    def compare_with_networkx(
        self,
        source: str,
        weight_key: str = "w_efficacy",
        targets: list[str] | None = None,
    ) -> dict:
        """
        Run both NetworkX Dijkstra and Julia DMY, compare results.

        Returns dict with distances from both methods and max discrepancy.
        """
        # NetworkX Dijkstra
        nx_dists = {}
        try:
            nx_all = nx.single_source_dijkstra_path_length(
                self.G, source, weight=weight_key
            )
            nx_dists = dict(nx_all)
        except nx.NodeNotFound:
            pass

        # Julia DMY
        dmy_result = self.run_sssp(source, weight_key, targets)

        # Compare
        comparison = []
        check_nodes = targets or list(nx_dists.keys())
        max_diff = 0.0

        for node_id in check_nodes:
            nx_d = nx_dists.get(node_id, float("inf"))
            dmy_d = dmy_result.distances.get(node_id, float("inf"))
            diff = abs(nx_d - dmy_d) if math.isfinite(nx_d) and math.isfinite(dmy_d) else 0.0
            max_diff = max(max_diff, diff)
            comparison.append({
                "node": node_id,
                "name": self.G.nodes[node_id].get("name", node_id),
                "networkx": nx_d,
                "dmy": dmy_d,
                "diff": diff,
                "match": diff < 1e-6,
            })

        return {
            "source": source,
            "weight_key": weight_key,
            "n_compared": len(comparison),
            "max_discrepancy": max_diff,
            "all_match": max_diff < 1e-6,
            "details": comparison,
            "dmy_result": dmy_result,
        }
