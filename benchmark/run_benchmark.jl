# run_benchmark.jl --- Experimental Evaluation harness (Section 7, requirement 1: Scale)
#
# Measures single-source shortest-path wall-clock time for three implementations
# on random sparse directed graphs (m approx 2n) at increasing scale, and reports
# the crossover point and maximum speedup.
#
# Three timed implementations (requirement 6, Baselines):
#   1. DMY               -- OptimShortestPaths.dmy_sssp!         (the algorithm under test)
#   2. Dijkstra (tuned)  -- internal hand-rolled binary-heap     (fair internal baseline)
#   3. Dijkstra (Graphs) -- Graphs.jl dijkstra_shortest_paths    (external sanity anchor)
#
# Methodology (requirement 5): BenchmarkTools.jl with warm-up, >=30 samples,
# median + bootstrap 95% CI, fixed RNG seed, machine/version stamped into the CSV
# header. GC time is recorded as a fraction so it can be accounted for.
#
# Usage:
#   julia --project=benchmark benchmark/run_benchmark.jl            # quick ladder (-> 1e4)
#   julia --project=benchmark benchmark/run_benchmark.jl --full     # full ladder (-> 1e5)
#   julia --project=benchmark benchmark/run_benchmark.jl --full --huge   # adds 3e5
#
# Output: benchmark/results/scale.csv

using OptimShortestPaths
using Graphs
using BenchmarkTools
using SparseArrays
using Random
using Statistics
using Printf
using InteractiveUtils: versioninfo

const OSP = OptimShortestPaths

# ---------------------------------------------------------------------------
# Tuned binary min-heap Dijkstra (lazy deletion; no per-decrease-key cost).
# Operates directly on the DMYGraph adjacency so it sees the exact same graph
# the DMY algorithm does -- a fair internal baseline.
# ---------------------------------------------------------------------------
mutable struct MinHeap
    data::Vector{Tuple{Float64,Int}}
    size::Int
end
MinHeap() = MinHeap(Tuple{Float64,Int}[], 0)

@inline function heap_push!(h::MinHeap, item::Tuple{Float64,Int})
    h.size += 1
    if h.size > length(h.data)
        push!(h.data, item)
    else
        @inbounds h.data[h.size] = item
    end
    i = h.size
    @inbounds while i > 1
        p = i >> 1
        if h.data[i][1] < h.data[p][1]
            h.data[i], h.data[p] = h.data[p], h.data[i]
            i = p
        else
            break
        end
    end
end

@inline function heap_pop!(h::MinHeap)
    @inbounds item = h.data[1]
    @inbounds h.data[1] = h.data[h.size]
    h.size -= 1
    i = 1
    @inbounds while true
        l = 2i; r = 2i + 1; smallest = i
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

function dijkstra_tuned(g::DMYGraph, src::Int)
    n = g.n_vertices
    dist = fill(Inf, n)
    dist[src] = 0.0
    done = falses(n)
    pq = MinHeap()
    heap_push!(pq, (0.0, src))
    while pq.size > 0
        d, u = heap_pop!(pq)
        done[u] && continue
        done[u] = true
        @inbounds for ei in g.adjacency_list[u]
            e = g.edges[ei]
            nd = d + g.weights[ei]
            if nd < dist[e.target]
                dist[e.target] = nd
                heap_push!(pq, (nd, e.target))
            end
        end
    end
    return dist
end

# ---------------------------------------------------------------------------
# Random sparse directed graph: n vertices, ~m_target unique edges (no self-loops,
# no parallel edges so all three implementations see an identical graph), with
# weights drawn uniformly from (0, 1]. Deterministic given `seed`.
# ---------------------------------------------------------------------------
function random_sparse_digraph(n::Int, m_target::Int; seed::Int=42)
    rng = MersenneTwister(seed)
    seen = Set{Tuple{Int,Int}}()
    srcs = Int[]; dsts = Int[]; wts = Float64[]
    # First lay down a random spanning-ish path so vertex 1 reaches a large
    # component (keeps the SSSP from terminating trivially on tiny frontiers).
    perm = randperm(rng, n)
    for i in 1:(n - 1)
        u, v = perm[i], perm[i + 1]
        push!(seen, (u, v)); push!(srcs, u); push!(dsts, v); push!(wts, rand(rng))
    end
    while length(srcs) < m_target
        u = rand(rng, 1:n); v = rand(rng, 1:n)
        (u == v || (u, v) in seen) && continue
        push!(seen, (u, v)); push!(srcs, u); push!(dsts, v); push!(wts, rand(rng))
    end
    edges = [OSP.Edge(srcs[i], dsts[i], i) for i in eachindex(srcs)]
    dmy_graph = DMYGraph(n, edges, wts)

    # Mirror into Graphs.jl: a SimpleDiGraph + sparse weight matrix (m nonzeros,
    # never the n^2 dense matrix -- critical at n=1e5).
    sg = SimpleDiGraph(n)
    for i in eachindex(srcs)
        add_edge!(sg, srcs[i], dsts[i])
    end
    distmx = sparse(srcs, dsts, wts, n, n)
    return dmy_graph, sg, distmx
end

# ---------------------------------------------------------------------------
# Correctness: all three must agree on distances (Inf == Inf for unreachable).
# ---------------------------------------------------------------------------
function distances_agree(a::Vector{Float64}, b::Vector{Float64}; atol=1e-9, rtol=1e-7)
    length(a) == length(b) || return false
    @inbounds for i in eachindex(a)
        ai, bi = a[i], b[i]
        if isinf(ai) || isinf(bi)
            (isinf(ai) && isinf(bi)) || return false
        elseif !isapprox(ai, bi; atol=atol, rtol=rtol)
            return false
        end
    end
    return true
end

# ---------------------------------------------------------------------------
# Bootstrap 95% CI of the median over BenchmarkTools sample times (ns).
# Deterministic given `seed`. Returns (median, lo, hi) in milliseconds.
# ---------------------------------------------------------------------------
function median_ci_ms(times_ns::Vector{Float64}; n_boot::Int=2000, seed::Int=1)
    rng = MersenneTwister(seed)
    k = length(times_ns)
    boot = Vector{Float64}(undef, n_boot)
    @inbounds for b in 1:n_boot
        acc = Vector{Float64}(undef, k)
        for j in 1:k
            acc[j] = times_ns[rand(rng, 1:k)]
        end
        boot[b] = median(acc)
    end
    med = median(times_ns)
    lo = quantile(boot, 0.025)
    hi = quantile(boot, 0.975)
    return (med / 1e6, lo / 1e6, hi / 1e6)
end

# ---------------------------------------------------------------------------
# Benchmark one callable on a single source. Returns NamedTuple of stats (ms).
# ---------------------------------------------------------------------------
function bench(f, args...; samples::Int=30, seconds::Float64=20.0)
    trial = @benchmark $f($(args)...) samples = samples seconds = seconds evals = 1
    times = Float64.(trial.times)
    med, lo, hi = median_ci_ms(times)
    gc_frac = isempty(trial.gctimes) ? 0.0 : mean(trial.gctimes) / mean(trial.times)
    return (median_ms=med, ci_lo_ms=lo, ci_hi_ms=hi,
            min_ms=minimum(times) / 1e6, n_samples=length(times), gc_frac=gc_frac)
end

# ---------------------------------------------------------------------------
# Main sweep.
# ---------------------------------------------------------------------------
function main()
    full = "--full" in ARGS
    huge = "--huge" in ARGS
    ns = full ? [1_000, 3_000, 10_000, 30_000, 100_000] : [1_000, 3_000, 10_000]
    huge && push!(ns, 300_000)
    src = 1
    seed = 42

    resultsdir = joinpath(@__DIR__, "results")
    mkpath(resultsdir)
    csvpath = joinpath(resultsdir, "scale.csv")

    open(csvpath, "w") do io
        # Methodology stamp (requirement 5).
        println(io, "# OptimShortestPaths Scale benchmark (Section 7, req 1)")
        println(io, "# julia=$(VERSION)  threads=$(Threads.nthreads())")
        println(io, "# cpu=$(Sys.cpu_info()[1].model)  ncpu=$(Sys.CPU_THREADS)  mem_gb=$(round(Sys.total_memory()/2^30, digits=1))")
        println(io, "# graph=random_sparse_digraph(m~2n, weights~U(0,1], no parallel/self edges)  seed=$seed  source=$src")
        println(io, "# timing=BenchmarkTools >=30 samples, median + bootstrap 95% CI, evals=1, warm-up via @benchmark tuning")
        println(io, "n,m,dmy_ms,dmy_lo,dmy_hi,dij_ms,dij_lo,dij_hi,graphs_ms,graphs_lo,graphs_hi,speedup_vs_tuned,speedup_vs_graphs,dmy_gc_frac,correct")
        flush(io)

        for n in ns
            m = 2 * n
            print(stderr, "n=$n m=$m  generating... ")
            dmy_graph, sg, distmx = random_sparse_digraph(n, m; seed=seed)

            # Correctness first (requirement: verify before timing).
            d_dmy = dmy_sssp!(dmy_graph, src)
            d_dij = dijkstra_tuned(dmy_graph, src)
            d_grf = dijkstra_shortest_paths(sg, src, distmx).dists
            ok = distances_agree(d_dmy, d_dij) && distances_agree(d_dmy, d_grf)
            print(stderr, ok ? "verified. " : "MISMATCH! ")

            print(stderr, "timing DMY... ")
            r_dmy = bench(dmy_sssp!, dmy_graph, src)
            print(stderr, "Dijkstra... ")
            r_dij = bench(dijkstra_tuned, dmy_graph, src)
            print(stderr, "Graphs.jl... ")
            r_grf = bench((g, s, w) -> dijkstra_shortest_paths(g, s, w), sg, src, distmx)

            sp_tuned = r_dij.median_ms / r_dmy.median_ms
            sp_graphs = r_grf.median_ms / r_dmy.median_ms

            @printf(io, "%d,%d,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.4f,%.4f,%.4f,%s\n",
                n, m,
                r_dmy.median_ms, r_dmy.ci_lo_ms, r_dmy.ci_hi_ms,
                r_dij.median_ms, r_dij.ci_lo_ms, r_dij.ci_hi_ms,
                r_grf.median_ms, r_grf.ci_lo_ms, r_grf.ci_hi_ms,
                sp_tuned, sp_graphs, r_dmy.gc_frac, ok)
            flush(io)
            @printf(stderr, "done.  DMY=%.3fms  Dij=%.3fms  Graphs=%.3fms  speedup(tuned)=%.2fx (graphs)=%.2fx\n",
                r_dmy.median_ms, r_dij.median_ms, r_grf.median_ms, sp_tuned, sp_graphs)
        end
    end

    println(stderr, "\nWrote $csvpath")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
