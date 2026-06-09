# Minimal, self-contained reproduction of the DMY correctness bug.
# Uses ONLY the package's own functions (dmy_sssp!, simple_dijkstra,
# compare_with_dijkstra). No external deps.
#
#   julia --project=. docs/dmy_bug_repro.jl

using OptimShortestPaths
using Random
using Printf
const O = OptimShortestPaths

# Random sparse digraph, m = 2n. `backbone` controls the spanning structure:
#   :tree -> forward spanning tree, parent = rand(1:i-1)  (what the test suite uses)
#   :path -> random Hamiltonian path over a random permutation (general graph)
function gen(n, m; backbone=:path, seed=42)
    rng = MersenneTwister(seed); seen = Set{Tuple{Int,Int}}()
    S = Int[]; D = Int[]; W = Float64[]
    add(u, v) = (push!(seen, (u, v)); push!(S, u); push!(D, v); push!(W, rand(rng) * 5 + 0.1))
    if backbone == :path
        perm = randperm(rng, n)
        for i in 1:(n - 1); add(perm[i], perm[i + 1]); end
    else
        for i in 2:n; add(rand(rng, 1:i - 1), i); end
    end
    while length(S) < m
        u = rand(rng, 1:n); v = rand(rng, 1:n)
        (u == v || (u, v) in seen) && continue
        add(u, v)
    end
    O.DMYGraph(n, [O.Edge(S[i], D[i], i) for i in eachindex(S)], W)
end

# Pure DMY recursive core, bypassing the Bellman-Ford safety net in dmy_sssp!.
function core_only(g, src)
    n = g.n_vertices; dist = fill(O.INF, n); parent = fill(0, n); dist[src] = 0.0
    O.recursive_layer!(g, dist, parent, collect(1:n), Set([src]), O.INF)
    dist
end

mism(a, b) = count(i -> !((isinf(a[i]) && isinf(b[i])) ||
                          (isfinite(a[i]) && isfinite(b[i]) && isapprox(a[i], b[i]; atol=1e-9))),
                   eachindex(a))

println("DMY correctness vs the package's own simple_dijkstra (ground truth), n=1000, m=2000\n")
println("backbone | compare_with_dijkstra | shipped dmy_sssp! errors | pure-core errors")
println("---------|-----------------------|--------------------------|-----------------")
for backbone in (:tree, :path)
    g = gen(1000, 2000; backbone=backbone, seed=42); src = 1
    truth = simple_dijkstra(g, src)
    shipped = dmy_sssp!(g, src)
    core = core_only(g, src)
    c = compare_with_dijkstra(g, src)   # the repo's own correctness comparator
    @printf("%-8s | match=%-5s disc=%-5d | %24d | %16d\n",
        backbone, string(c["results_match"]), length(c["discrepancies"]),
        mism(shipped, truth), mism(core, truth))
end

println("\nExpected: tree -> correct (what tests cover); path -> WRONG (general graphs).")
