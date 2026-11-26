using TDARipserer, Test, Random, Distances, SparseArrays

endpoints(diagrams) = [sort!([(birth(i), death(i)) for i in d]) for d in diagrams]

@testset "Exact reduced Rips barcodes" begin
    rng = MersenneTwister(230716333)
    grid = [(x, y, z) for x in 0:2 for y in 0:2 for z in 0:2]
    fixtures = [
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [Tuple(rand(rng, 2)) for _ in 1:45],
        [Tuple(rand(rng, 3)) for _ in 1:60],
        grid,
        shuffle(rng, grid),
        [(i, 0, 0) for i in 1:40],
    ]
    for points in fixtures
        complete = ripserer(points; alg=:reduced_rips)
        # The existing reducer on these same pointwise grades is the exact oracle.
        weights = Matrix(TDARipserer.distance_matrix(complete[1].filtration))
        for limit in (nothing, 0.0, 0.4, 1.0, sqrt(2.0), 3.0, Inf), cutoff in (0, 0.1)
            actual = ripserer(points; alg=:reduced_rips, threshold=limit, cutoff)
            oracle = ripserer(weights; threshold=limit, cutoff)
            @test endpoints(actual) == endpoints(oracle)
            @test threshold.(actual) == threshold.(oracle)
        end
    end

    square = first(fixtures)
    @test endpoints(ripserer(square; alg=:reduced_rips))[2] == [(1.0, sqrt(2.0))]
    @test endpoints(ripserer(square; alg=:reduced_rips, threshold=prevfloat(sqrt(2.0))))[2] ==
        [(1.0, Inf)]
    @test endpoints(ripserer(square; alg=:reduced_rips, threshold=1.0))[2] == [(1.0, Inf)]
    @test isempty(ripserer(square; alg=:reduced_rips, threshold=prevfloat(1.0))[2])
    @test all(isempty, ripserer(square; alg=:reduced_rips, cutoff=Inf))
    @test endpoints(ripserer(square; alg=:reduced_rips, dim_max=0)) ==
        endpoints(ripserer(square; dim_max=0))
    typed = ripserer(Rips{Int32}, square; alg=:reduced_rips)
    interval = only(typed[2])
    @test birth_simplex(interval) isa Simplex{1,Float64,Int32}
    @test death_simplex(interval) isa Simplex{2,Float64,Int32}
    @test birth(death_simplex(interval)) == death(interval)

    shifted = [(3.0 + x, y) for (x, y) in square]
    @test endpoints(ripserer(vcat(square, shifted); alg=:reduced_rips, threshold=1.0))[2] ==
        [(1.0, Inf), (1.0, Inf)]
    @test_logs (:warn, "points not unique") ripserer([square; [(-0.0, 0.0)]];
        alg=:reduced_rips)
    @test endpoints(ripserer([(0f0, 0f0), (1f0, 0f0)]; alg=:reduced_rips)) ==
        [[(0.0, 1.0), (0.0, Inf)], []]
    @test all(isempty, ripserer(NTuple{3,Float64}[]; alg=:reduced_rips))
    @test endpoints(ripserer([(0.0, 0.0)]; alg=:reduced_rips)) == [[(0.0, Inf)], []]

    @test_throws ArgumentError ripserer(square; alg=:reduced_rips, dim_max=2)
    @test_throws ArgumentError ripserer(square; alg=:reduced_rips, modulus=3)
    @test_throws ArgumentError ripserer(square; alg=:reduced_rips, reps=true)
    @test_throws ArgumentError ripserer(square; alg=:reduced_rips, cutoff=-1)
    @test_throws ArgumentError ripserer(square; alg=:reduced_rips, metric=Cityblock())
    @test_throws ArgumentError ripserer(square; alg=:reduced_rips, threshold=NaN)
    @test_throws ArgumentError ripserer([(NaN, 0.0)]; alg=:reduced_rips)
    @test_throws ArgumentError ripserer([(0, 0, 0, 0)]; alg=:reduced_rips)
    @test_throws ArgumentError ripserer(Rips(square); alg=:reduced_rips)
    @test_throws ArgumentError ripserer(TDARipserer.distances(square); alg=:reduced_rips)
end

@testset "Direct radius graph" begin
    for T in (Float32, Float64)
        points = NTuple{3,T}[(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (4, 0, 0)]
        for limit in (zero(T), prevfloat(one(T)), one(T), sqrt(T(2)))
            direct = Rips(points; threshold=limit, radius_graph=true)
            @test issparse(TDARipserer.adjacency_matrix(direct))
            baseline = Rips(points; threshold=limit, sparse=true)
            @test TDARipserer.adjacency_matrix(direct) == TDARipserer.adjacency_matrix(baseline)
            @test endpoints(ripserer(direct; dim_max=2)) ==
                endpoints(ripserer(baseline; dim_max=2))
        end
        @test endpoints(ripserer(points; radius_graph=true, threshold=1, metric=Cityblock())) ==
            endpoints(ripserer(points; sparse=true, threshold=1, metric=Cityblock()))
    end
end
