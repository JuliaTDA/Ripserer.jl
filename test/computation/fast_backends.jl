using TDARipserer
using Distances: Cityblock, Euclidean, evaluate
using Random
using SparseArrays
using Test

const PH = TDARipserer
_signature(diagrams) = [sort!([(birth(b), death(b)) for b in d]) for d in diagrams]

@testset "Edge collapse with short thresholds and vertex grades" begin
    edgeless = Rips([(0.0,), (1.0,)]; radius_graph=true, threshold=0.5)
    @test _signature(ripserer(EdgeCollapsedRips(edgeless); dim_max=0)) ==
        [[(0.0, Inf), (0.0, Inf)]]
    for sparse_input in (false, true)
        partial = Rips([1.0 4.0; 4.0 3.0]; threshold=2.0, sparse=sparse_input)
        @test _signature(ripserer(EdgeCollapsedRips(partial); dim_max=0)) == [[(1.0, Inf)]]
    end
    for D in ([1.0 4.0; 4.0 3.0],
              [-4.0 -2 -1; -2 -3 -2; -1 -2 -2],
              [-2.0 0.5 1.5 2.5; 0.5 -1 2 2; 1.5 2 0 1; 2.5 2 1 1])
        for threshold in (-3.0, -2.0, -1.0, 0.25, 1.5, 2.5, 4.0), adj in (D, sparse(D))
            f = Rips(adj; threshold)
            collapsed = EdgeCollapsedRips(f)
            @test PH.births(collapsed) == PH.births(f)
            @test EdgeCollapsedRips(collapsed).adj === collapsed.adj
            @test PH.threshold(collapsed) == PH.threshold(f)
            dim_max = size(D, 1) - 2
            @test _signature(ripserer(f; dim_max)) == _signature(ripserer(collapsed; dim_max))
        end
    end
end

@testset "Indexed radius graphs" begin
    for T in (Float32, Float64)
        points = [(T(i), T(j)) for i in 1:4 for j in 1:4]
        for radius in T[0, prevfloat(one(T)), 1, sqrt(T(2)), 10]
            graph = Rips(points; radius_graph=true, threshold=radius)
            @test issparse(graph.adj)
            metric = Euclidean(1e-12)
            oracle = [let d = evaluate(metric, points[i], points[j])
                        d <= radius ? d : zero(T)
                      end for i in eachindex(points), j in eachindex(points)]
            @test Matrix(graph.adj) == oracle
            @test _signature(ripserer(points; radius_graph=true, threshold=radius, dim_max=2)) ==
                _signature(ripserer(oracle |> sparse; threshold=radius, dim_max=2))
        end
    end
    points = [(0.3, 0.5), (0.6, 0.1), (1.0, 2.0)]
    r = evaluate(Euclidean(1e-12), points[1], points[2])
    @test Rips(points; threshold=r, radius_graph=true).adj[1, 2] == r
    @test iszero(Rips(points; threshold=prevfloat(r), radius_graph=true).adj[1, 2])
    @test Rips(points; radius_graph=true, threshold=1, metric=Cityblock()).adj ==
        Rips(points; sparse=true, threshold=1, metric=Cityblock()).adj
end

@testset "RedZeD endpoint and simplex parity" begin
    rng = MersenneTwister(19)
    square = [0.0 1 2 1; 1 0 1 2; 2 1 0 1; 1 2 1 0]
    @test _signature(ripserer(Rips(square); dim_max=9, alg=:redzed))[2] == [(1, 2)]
    @test _signature(ripserer(Rips(zeros(1, 1)); dim_max=0, alg=:redzed)) == [[(0, Inf)]]
    for T in (Float32, Float64), fixture in (:ties, :metric, :nonmetric)
        D = if fixture == :ties
            T.(square)
        elseif fixture == :metric
            PH.distances([Tuple(p) for p in eachcol(rand(rng, T, 4, 24))])
        else
            A = rand(rng, T, 24, 24)
            [i == j ? zero(T) : A[min(i, j), max(i, j)] for i in 1:24, j in 1:24]
        end
        for limit in T[0, 0.7, 1, Inf], cutoff in (0.0, 0.3)
            dim = min(3, size(D, 1) - 2)
            f = Rips(D; threshold=limit)
            actual = ripserer(f; alg=:redzed, dim_max=dim, cutoff)
            @test _signature(actual) == _signature(ripserer(f; dim_max=dim, cutoff))
            for bar in Iterators.flatten(actual)
                @test birth(birth_simplex(bar)) == birth(bar)
                @test death_simplex(bar) === nothing ? !isfinite(bar) :
                    birth(death_simplex(bar)) == death(bar)
            end
        end
    end
    @test_throws ArgumentError ripserer(square; alg=:redzed, reps=true)
    @test_throws ArgumentError ripserer(square; alg=:redzed, modulus=3)
    @test_throws ArgumentError ripserer(sparse(square); alg=:redzed)
    @test_throws ArgumentError ripserer(square; alg=:redzed, cutoff=-1)
end

@testset "Cubical union-find with tied grades and truncated grids" begin
    rng = MersenneTwister(47)
    for shape in ((2, 3), (4, 5), (1, 8), (3, 4, 5), (3, 4, 1))
        for _ in 1:6, t in (-3, -1, 0, 2)
            f = Cubical(rand(rng, -2:2, shape...); threshold=t)
            dim = length(shape) - 1
            actual = ripserer(f; alg=:cubical_union_find, dim_max=dim)
            @test _signature(actual) == _signature(ripserer(f; dim_max=dim))
            for bar in Iterators.flatten(actual)
                @test birth(birth_simplex(bar)) == birth(bar)
                @test death_simplex(bar) === nothing ? !isfinite(bar) :
                    birth(death_simplex(bar)) == death(bar)
            end
        end
    end
    hollow = zeros(Int, 3, 3, 3)
    hollow[2, 2, 2] = 1
    @test _signature(ripserer(Cubical, hollow; alg=:cubical_union_find, dim_max=2))[3] == [(0, 1)]
    @test _signature(ripserer(Cubical, hollow; alg=:cubical_union_find, dim_max=2, threshold=0))[3] ==
        [(0, Inf)]
    @test _signature(ripserer(Cubical, hollow; alg=:cubical_union_find, reps=true)) ==
        _signature(ripserer(Cubical, hollow; reps=true))
end

@testset "Automatic radius graph selection" begin
    for T in (Float32, Float64)
        points = [(T(i), T(i % 2)) for i in 1:512]
        for radius in (prevfloat(sqrt(T(2))), sqrt(T(2)), nextfloat(sqrt(T(2))))
            indexed = Rips(points; radius_graph=true, threshold=radius)
            automatic = Rips(points; radius_graph=:auto, threshold=radius)
            pairwise = Rips(points; radius_graph=false, sparse=true, threshold=radius)
            @test issparse(automatic.adj)
            @test automatic.adj == indexed.adj == pairwise.adj
            @test _signature(ripserer(automatic)) == _signature(ripserer(indexed))
        end
        @test Rips(points; radius_graph=:auto, threshold=T(1024)).adj ==
            Rips(points; radius_graph=false, threshold=T(1024)).adj
        @test !issparse(Rips(points; radius_graph=:auto, threshold=T(1024)).adj)
        @test Rips(points; metric=Cityblock(), radius_graph=:auto, threshold=1).adj ==
            Rips(points; metric=Cityblock(), threshold=1).adj
    end
    @test_throws ArgumentError Rips([(0.0,), (1.0,)]; radius_graph=:unknown)
end

@testset "Parallel apparent pairs" begin
    rng = MersenneTwister(812)
    for n in (8, 12, 24), format in (:dense, :sparse, :collapse), threshold in (2, 3, 5)
        distances = zeros(Int, n, n)
        for j in 2:n, i in 1:(j - 1)
            distances[i, j] = distances[j, i] = rand(rng, 1:5)
        end
        f = Rips(distances; threshold, sparse=format == :sparse)
        f = format == :collapse ? EdgeCollapsedRips(f) : f
        @test _signature(ripserer(f; dim_max=3)) ==
            _signature(ripserer(f; dim_max=3, alg=:parallel_cohomology))
    end
    @test_throws ArgumentError ripserer([(0.0,), (1.0,)]; alg=:parallel_cohomology, modulus=3)
end
