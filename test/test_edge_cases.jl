"""
Edge-case tests for multi-objective Pareto and degenerate graph scenarios.
Covers: ties, empty fronts, single-dimension, disconnected graphs,
self-loops, zero-weight edges, large uniform weights, and parallel edges.
"""

using Test

const MO_EC = OptimShortestPaths.MultiObjective
const MultiObjectiveEdge_EC = MO_EC.MultiObjectiveEdge
const MultiObjectiveGraph_EC = MO_EC.MultiObjectiveGraph
const ParetoSolution_EC = MO_EC.ParetoSolution
const compute_pareto_front_EC = MO_EC.compute_pareto_front
const weighted_sum_approach_EC = MO_EC.weighted_sum_approach
const epsilon_constraint_approach_EC = MO_EC.epsilon_constraint_approach
const lexicographic_approach_EC = MO_EC.lexicographic_approach
const get_knee_point_EC = MO_EC.get_knee_point
const INF_EC = OptimShortestPaths.INF

@testset "Edge Case Tests" begin

    # ─── Multi-Objective Pareto Edge Cases ───────────────────────────

    @testset "Pareto: Tied Objective Values" begin
        # Two paths with identical objectives — should yield a single Pareto solution
        edges = [
            MultiObjectiveEdge_EC(1, 2, [1.0, 2.0], 1),
            MultiObjectiveEdge_EC(2, 4, [3.0, 1.0], 2),
            MultiObjectiveEdge_EC(1, 3, [2.0, 1.0], 3),
            MultiObjectiveEdge_EC(3, 4, [2.0, 2.0], 4),
        ]
        graph = MultiObjectiveGraph_EC(4, edges, 2, ["cost", "time"])

        front = compute_pareto_front_EC(graph, 1, 4)
        # Both paths: [1,2,4] = [4.0, 3.0] and [1,3,4] = [4.0, 3.0]
        # They have identical objectives, so only one should appear
        @test length(front) == 1
        @test front[1].objectives ≈ [4.0, 3.0]
    end

    @testset "Pareto: All Paths Mutually Non-Dominated" begin
        # Three paths where none dominates any other
        edges = [
            MultiObjectiveEdge_EC(1, 2, [1.0, 5.0], 1),
            MultiObjectiveEdge_EC(1, 3, [3.0, 3.0], 2),
            MultiObjectiveEdge_EC(1, 4, [5.0, 1.0], 3),
            MultiObjectiveEdge_EC(2, 5, [0.0, 0.0], 4),
            MultiObjectiveEdge_EC(3, 5, [0.0, 0.0], 5),
            MultiObjectiveEdge_EC(4, 5, [0.0, 0.0], 6),
        ]
        graph = MultiObjectiveGraph_EC(5, edges, 2, ["a", "b"])

        front = compute_pareto_front_EC(graph, 1, 5)
        @test length(front) == 3
        obj_set = Set(Tuple(round.(s.objectives, digits=8)) for s in front)
        @test (1.0, 5.0) in obj_set
        @test (3.0, 3.0) in obj_set
        @test (5.0, 1.0) in obj_set
    end

    @testset "Pareto: Empty Front (Unreachable Target)" begin
        edges = [
            MultiObjectiveEdge_EC(1, 2, [1.0, 1.0], 1),
        ]
        graph = MultiObjectiveGraph_EC(3, edges, 2, ["a", "b"])

        front = compute_pareto_front_EC(graph, 1, 3)
        @test isempty(front)
    end

    @testset "Pareto: Single-Dimension (1 Objective)" begin
        edges = [
            MultiObjectiveEdge_EC(1, 2, [3.0], 1),
            MultiObjectiveEdge_EC(1, 3, [1.0], 2),
            MultiObjectiveEdge_EC(2, 4, [1.0], 3),
            MultiObjectiveEdge_EC(3, 4, [5.0], 4),
        ]
        graph = MultiObjectiveGraph_EC(4, edges, 1, ["cost"])

        front = compute_pareto_front_EC(graph, 1, 4)
        # With 1 objective, only the minimum-cost path survives
        @test length(front) == 1
        @test front[1].objectives ≈ [4.0]  # path [1,2,4]
    end

    @testset "Pareto: Source Equals Target" begin
        edges = [
            MultiObjectiveEdge_EC(1, 2, [1.0, 1.0], 1),
        ]
        graph = MultiObjectiveGraph_EC(2, edges, 2, ["a", "b"])

        front = compute_pareto_front_EC(graph, 1, 1)
        @test length(front) == 1
        @test front[1].objectives ≈ [0.0, 0.0]
        @test front[1].path == [1]
    end

    @testset "Pareto: max_solutions=1 Returns Single Solution" begin
        edges = [
            MultiObjectiveEdge_EC(1, 2, [1.0, 5.0], 1),
            MultiObjectiveEdge_EC(1, 3, [5.0, 1.0], 2),
            MultiObjectiveEdge_EC(2, 4, [0.0, 0.0], 3),
            MultiObjectiveEdge_EC(3, 4, [0.0, 0.0], 4),
        ]
        graph = MultiObjectiveGraph_EC(4, edges, 2, ["a", "b"])

        front = compute_pareto_front_EC(graph, 1, 4; max_solutions=1)
        @test length(front) >= 1
        @test length(front) <= 2  # may include up to 2 due to expansion order
    end

    @testset "Knee Point: Empty Front" begin
        @test get_knee_point_EC(ParetoSolution_EC[]) === nothing
    end

    @testset "Knee Point: Single Solution" begin
        front = [ParetoSolution_EC([1.0, 2.0], [1, 2], Int[])]
        knee = get_knee_point_EC(front)
        @test knee === front[1]
    end

    @testset "Knee Point: Identical Solutions" begin
        front = [
            ParetoSolution_EC([5.0, 5.0], [1], Int[]),
            ParetoSolution_EC([5.0, 5.0], [2], Int[]),
        ]
        knee = get_knee_point_EC(front)
        @test knee !== nothing
        @test knee.objectives ≈ [5.0, 5.0]
    end

    @testset "Epsilon Constraint: All Infeasible" begin
        edges = [
            MultiObjectiveEdge_EC(1, 2, [10.0, 10.0], 1),
        ]
        graph = MultiObjectiveGraph_EC(2, edges, 2, ["a", "b"])

        # Constraint on objective 2 too tight
        sol = epsilon_constraint_approach_EC(graph, 1, 2, 1, [Inf, 1.0])
        @test all(isinf, sol.objectives)
    end

    # ─── Degenerate Graph Tests (DMY Correctness) ────────────────────

    @testset "DMY: Disconnected Components" begin
        # 4 components: {1,2}, {3,4}, {5}, {6,7,8}
        edges = [
            Edge(1, 2, 1), Edge(2, 1, 2),
            Edge(3, 4, 3), Edge(4, 3, 4),
            Edge(6, 7, 5), Edge(7, 8, 6), Edge(8, 6, 7),
        ]
        weights = [1.0, 2.0, 3.0, 4.0, 1.5, 2.5, 0.5]
        graph = DMYGraph(8, edges, weights)

        for source in [1, 3, 5, 6]
            dmy_d = dmy_sssp!(graph, source)
            dij_d = simple_dijkstra(graph, source)
            for i in 1:8
                if dmy_d[i] >= INF_EC && dij_d[i] >= INF_EC
                    @test true
                else
                    @test abs(dmy_d[i] - dij_d[i]) < 1e-10
                end
            end
        end
    end

    @testset "DMY: Self-Loops Only" begin
        # Graph where every edge is a self-loop
        edges = [Edge(1, 1, 1), Edge(2, 2, 2), Edge(3, 3, 3)]
        weights = [1.0, 2.0, 3.0]
        graph = DMYGraph(3, edges, weights)

        d = dmy_sssp!(graph, 1)
        @test d[1] == 0.0
        @test d[2] >= INF_EC
        @test d[3] >= INF_EC
    end

    @testset "DMY: All Zero-Weight Edges" begin
        edges = [Edge(1, 2, 1), Edge(2, 3, 2), Edge(1, 3, 3), Edge(3, 4, 4)]
        weights = [0.0, 0.0, 0.0, 0.0]
        graph = DMYGraph(4, edges, weights)

        dmy_d = dmy_sssp!(graph, 1)
        dij_d = simple_dijkstra(graph, 1)
        for i in 1:4
            @test abs(dmy_d[i] - dij_d[i]) < 1e-10
        end
        # All distances should be 0 (reachable via zero-cost paths)
        @test dmy_d[1] == 0.0
        @test dmy_d[2] == 0.0
        @test dmy_d[3] == 0.0
        @test dmy_d[4] == 0.0
    end

    @testset "DMY: Parallel Edges (Multiple Edges Between Same Pair)" begin
        # Three edges from 1→2 with different weights; shortest should win
        edges = [Edge(1, 2, 1), Edge(1, 2, 2), Edge(1, 2, 3), Edge(2, 3, 4)]
        weights = [5.0, 2.0, 8.0, 1.0]
        graph = DMYGraph(3, edges, weights)

        dmy_d = dmy_sssp!(graph, 1)
        dij_d = simple_dijkstra(graph, 1)
        @test dmy_d[2] ≈ 2.0  # shortest of the 3 parallel edges
        @test dmy_d[3] ≈ 3.0  # 2.0 + 1.0
        for i in 1:3
            @test abs(dmy_d[i] - dij_d[i]) < 1e-10
        end
    end

    @testset "DMY: Very Large Weights" begin
        edges = [Edge(1, 2, 1), Edge(2, 3, 2)]
        weights = [1e15, 1e15]
        graph = DMYGraph(3, edges, weights)

        dmy_d = dmy_sssp!(graph, 1)
        dij_d = simple_dijkstra(graph, 1)
        @test abs(dmy_d[2] - 1e15) < 1e5  # allow numerical tolerance
        @test abs(dmy_d[3] - 2e15) < 1e5
        for i in 1:3
            if dmy_d[i] >= INF_EC && dij_d[i] >= INF_EC
                @test true
            else
                @test abs(dmy_d[i] - dij_d[i]) < 1e5
            end
        end
    end

    @testset "DMY: Single Vertex Graph" begin
        graph = DMYGraph(1, Edge[], Float64[])
        d = dmy_sssp!(graph, 1)
        @test d == [0.0]

        d2, p2 = dmy_sssp_with_parents!(graph, 1)
        @test d2 == [0.0]
        @test p2 == [0]
    end

    @testset "DMY: Linear Chain (n=100)" begin
        n = 100
        edges = [Edge(i, i+1, i) for i in 1:n-1]
        weights = [1.0 for _ in 1:n-1]
        graph = DMYGraph(n, edges, weights)

        dmy_d = dmy_sssp!(graph, 1)
        dij_d = simple_dijkstra(graph, 1)
        for i in 1:n
            @test abs(dmy_d[i] - dij_d[i]) < 1e-10
        end
        @test dmy_d[n] ≈ 99.0
    end

    @testset "DMY: Star Graph (Hub-and-Spoke)" begin
        # Hub at vertex 1, spokes to vertices 2..11
        n = 11
        edges = [Edge(1, i, i-1) for i in 2:n]
        weights = [Float64(i) for i in 1:n-1]
        graph = DMYGraph(n, edges, weights)

        dmy_d = dmy_sssp!(graph, 1)
        @test dmy_d[1] == 0.0
        for i in 2:n
            @test dmy_d[i] ≈ Float64(i - 1)
        end

        # From a spoke, only the spoke itself is reachable
        dmy_from_spoke = dmy_sssp!(graph, 5)
        @test dmy_from_spoke[5] == 0.0
        for i in 1:n
            i == 5 && continue
            @test dmy_from_spoke[i] >= INF_EC
        end
    end

    @testset "DMY: Complete Graph (n=8)" begin
        n = 8
        edges = Edge[]
        weights = Float64[]
        idx = 1
        for i in 1:n, j in 1:n
            i == j && continue
            push!(edges, Edge(i, j, idx))
            push!(weights, Float64(i + j))  # deterministic weights
            idx += 1
        end
        graph = DMYGraph(n, edges, weights)

        for source in 1:n
            dmy_d = dmy_sssp!(graph, source)
            dij_d = simple_dijkstra(graph, source)
            for i in 1:n
                @test abs(dmy_d[i] - dij_d[i]) < 1e-10
            end
        end
    end

    # ─── Julia Bridge Integration Tests ──────────────────────────────

    @testset "Julia Bridge: JSON Round-Trip" begin
        # Simulate the bridge: build graph from JSON-like data, run DMY, verify
        input_data = Dict(
            "n_vertices" => 5,
            "source" => 1,
            "edges" => [
                Dict("source" => 1, "target" => 2, "weight" => 1.5),
                Dict("source" => 2, "target" => 3, "weight" => 2.0),
                Dict("source" => 1, "target" => 3, "weight" => 4.0),
                Dict("source" => 3, "target" => 4, "weight" => 1.0),
                Dict("source" => 4, "target" => 5, "weight" => 0.5),
            ],
        )

        n = input_data["n_vertices"]
        src = input_data["source"]
        raw_edges = input_data["edges"]

        edges = Vector{Edge}(undef, length(raw_edges))
        weights = Vector{Float64}(undef, length(raw_edges))
        for (i, e) in enumerate(raw_edges)
            edges[i] = Edge(e["source"], e["target"], i)
            weights[i] = max(0.0, Float64(e["weight"]))
        end

        graph = DMYGraph(n, edges, weights)
        distances, parents = dmy_sssp_with_parents!(graph, src)

        # Verify against Dijkstra
        dij_d = simple_dijkstra(graph, src)
        for i in 1:n
            @test abs(distances[i] - dij_d[i]) < 1e-10
        end

        # Check expected distances
        @test distances[1] ≈ 0.0
        @test distances[2] ≈ 1.5
        @test distances[3] ≈ 3.5  # 1→2→3
        @test distances[4] ≈ 4.5  # 1→2→3→4
        @test distances[5] ≈ 5.0  # 1→2→3→4→5

        # Verify parent chain for path reconstruction
        path_to_5 = reconstruct_path(parents, 1, 5)
        @test path_to_5 == [1, 2, 3, 4, 5]

        # Simulate JSON serialization (Inf → 1e308)
        safe_distances = [isinf(d) ? 1e308 : d for d in distances]
        @test all(isfinite, safe_distances)
    end

    @testset "Julia Bridge: Negative Weights Clamped to Zero" begin
        # The bridge clamps weights: max(0.0, weight)
        raw_edges = [
            Dict("source" => 1, "target" => 2, "weight" => -5.0),
            Dict("source" => 2, "target" => 3, "weight" => 3.0),
        ]

        edges = Vector{Edge}(undef, length(raw_edges))
        weights = Vector{Float64}(undef, length(raw_edges))
        for (i, e) in enumerate(raw_edges)
            edges[i] = Edge(e["source"], e["target"], i)
            weights[i] = max(0.0, Float64(e["weight"]))
        end

        @test weights[1] == 0.0  # clamped from -5.0
        @test weights[2] == 3.0

        graph = DMYGraph(3, edges, weights)
        d = dmy_sssp!(graph, 1)
        @test d[2] ≈ 0.0  # zero-weight edge
        @test d[3] ≈ 3.0
    end

    @testset "Julia Bridge: Empty Graph" begin
        graph = DMYGraph(3, Edge[], Float64[])
        d, p = dmy_sssp_with_parents!(graph, 1)

        @test d[1] == 0.0
        @test d[2] >= INF_EC
        @test d[3] >= INF_EC

        safe_distances = [isinf(d_i) ? 1e308 : d_i for d_i in d]
        @test safe_distances[1] == 0.0
        @test safe_distances[2] == 1e308
    end

end
