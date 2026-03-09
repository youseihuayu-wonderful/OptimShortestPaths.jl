
using OptimShortestPaths
using JSON

# ---------- Binary min-heap priority queue ----------
mutable struct MinHeap
    data::Vector{Tuple{Float64,Int}}
    size::Int
end
MinHeap() = MinHeap(Tuple{Float64,Int}[], 0)

function heap_push!(h::MinHeap, item::Tuple{Float64,Int})
    h.size += 1
    if h.size > length(h.data)
        push!(h.data, item)
    else
        h.data[h.size] = item
    end
    # Sift up
    i = h.size
    while i > 1
        p = i >> 1
        if h.data[i][1] < h.data[p][1]
            h.data[i], h.data[p] = h.data[p], h.data[i]
            i = p
        else
            break
        end
    end
end

function heap_pop!(h::MinHeap)
    item = h.data[1]
    h.data[1] = h.data[h.size]
    h.size -= 1
    # Sift down
    i = 1
    while true
        l = 2i
        r = 2i + 1
        smallest = i
        if l <= h.size && h.data[l][1] < h.data[smallest][1]
            smallest = l
        end
        if r <= h.size && h.data[r][1] < h.data[smallest][1]
            smallest = r
        end
        if smallest != i
            h.data[i], h.data[smallest] = h.data[smallest], h.data[i]
            i = smallest
        else
            break
        end
    end
    return item
end

# ---------- Proper Dijkstra with binary heap ----------
function dijkstra_heap!(adj_list, edges_arr, weights_arr, n, src)
    dist = fill(Inf, n)
    dist[src] = 0.0
    visited = falses(n)
    pq = MinHeap()
    heap_push!(pq, (0.0, src))
    while pq.size > 0
        d, u = heap_pop!(pq)
        visited[u] && continue
        visited[u] = true
        dist[u] = d
        for edge_idx in adj_list[u]
            e = edges_arr[edge_idx]
            w = weights_arr[edge_idx]
            nd = d + w
            if nd < dist[e.target]
                dist[e.target] = nd
                heap_push!(pq, (nd, e.target))
            end
        end
    end
    return dist
end

function main()
    input_path = ARGS[1]
    output_path = ARGS[2]

    data = JSON.parsefile(input_path)
    n = data["n_vertices"]
    edge_list = data["edges"]
    sources = data["sources"]
    n_warmup = get(data, "n_warmup", 1)

    # Build graph
    edges = Edge[]
    weights = Float64[]
    for e in edge_list
        idx = length(edges) + 1
        push!(edges, Edge(e["source"], e["target"], idx))
        push!(weights, e["weight"])
    end
    g = DMYGraph(n, edges, weights)

    # Warmup both algorithms (JIT compilation)
    for _ in 1:n_warmup
        dmy_sssp_with_parents!(g, sources[1])
        dijkstra_heap!(g.adjacency_list, g.edges, g.weights, n, sources[1])
    end

    # Timed DMY runs
    times_ns = Int64[]
    all_dists = Vector{Float64}[]
    for src in sources
        t0 = time_ns()
        dists, _ = dmy_sssp_with_parents!(g, src)
        t1 = time_ns()
        push!(times_ns, t1 - t0)
        push!(all_dists, dists)
    end

    # Timed Dijkstra runs (proper binary heap)
    dijkstra_times_ns = Int64[]
    for src in sources
        t0 = time_ns()
        dijkstra_heap!(g.adjacency_list, g.edges, g.weights, n, src)
        t1 = time_ns()
        push!(dijkstra_times_ns, t1 - t0)
    end

    # Convert Inf to large number for JSON
    safe_dists = [
        [d < 1e300 ? d : 1e308 for d in ds]
        for ds in all_dists
    ]

    result = Dict(
        "dmy_times_ns" => times_ns,
        "dijkstra_times_ns" => dijkstra_times_ns,
        "n_vertices" => n,
        "n_edges" => length(edges),
        "n_queries" => length(sources),
        "distances" => safe_dists,
    )

    open(output_path, "w") do f
        JSON.print(f, result)
    end
end

main()
