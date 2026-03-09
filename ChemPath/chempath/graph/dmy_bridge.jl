"""
Julia bridge script for ChemPath → OptimShortestPaths.jl integration.

Reads a JSON edge list from Python, builds a DMYGraph, runs dmy_sssp!,
and writes distances + parent pointers back to JSON.

Usage (called automatically by julia_bridge.py):
    julia --project=/path/to/OptimShortestPaths.jl dmy_bridge.jl input.json output.json
"""

using JSON
using OptimShortestPaths

function main()
    if length(ARGS) < 2
        error("Usage: julia dmy_bridge.jl <input.json> <output.json>")
    end

    input_path = ARGS[1]
    output_path = ARGS[2]

    # Read input
    data = JSON.parsefile(input_path)

    n_vertices = data["n_vertices"]
    source = data["source"]

    # Build edges and weights
    raw_edges = data["edges"]
    n_edges = length(raw_edges)

    edges = Vector{Edge}(undef, n_edges)
    weights = Vector{Float64}(undef, n_edges)

    for (i, e) in enumerate(raw_edges)
        edges[i] = Edge(e["source"], e["target"], i)
        weights[i] = max(0.0, Float64(e["weight"]))
    end

    # Build graph
    graph = DMYGraph(n_vertices, edges, weights)

    # Run DMY SSSP
    distances, parents = dmy_sssp_with_parents!(graph, source)

    # Replace Inf with a large finite number for JSON serialization
    safe_distances = [isinf(d) ? 1e308 : d for d in distances]

    # Write output
    result = Dict(
        "n_vertices" => n_vertices,
        "n_edges" => n_edges,
        "source" => source,
        "distances" => safe_distances,
        "parents" => parents,
        "algorithm" => "DMY O(m log^(2/3) n)",
    )

    open(output_path, "w") do f
        JSON.print(f, result)
    end
end

main()
