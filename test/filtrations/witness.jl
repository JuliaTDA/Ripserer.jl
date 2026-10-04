using TDARipserer
using TDAPersistenceDiagrams
using Test

@testset "Witness filtrations" begin
    # A square with witnesses at landmarks and edge midpoints has a loop
    # before diagonal edges/triangles enter.
    points = [[cos(t), sin(t)] for t in range(0, 2π; length=65)[1:64]]
    landmarks = [1, 17, 33, 49]
    for F in (Witness, LazyWitness)
        flt = F(points, landmarks)
        @test flt.landmarks == landmarks
        @test TDARipserer.nv(flt) == 4
        @test TDARipserer.births(flt) == zeros(4)
        diagrams = ripserer(flt)
        @test length(diagrams[2]) == 1
        @test persistence(only(diagrams[2])) > 0.1
        @test diagrams == ripserer(flt; alg=:homology)
        @test diagrams == ripserer(flt; alg=:involuted)
        @test length(representative(only(ripserer(flt; reps=true)[2]))) > 0
        @test ripserer(F, points, landmarks) == diagrams
        @test_throws ArgumentError F(points, [1, 1])
        @test_throws ArgumentError F(points, [0])
        @test_throws ArgumentError F(points, Int[])
        @test_throws ArgumentError F(zeros(0, 2))
        @test_throws ArgumentError F([-1.0 0.0])
        @test_throws ArgumentError F(points, landmarks; threshold=-1)
    end

    d = [0.0 1.0 2.0; 1.0 0.0 1.0; 2.0 1.0 0.0]
    lazy = LazyWitness(d; nu=0)
    @test TDARipserer.adjacency_matrix(lazy) == [0 1 1; 1 0 1; 1 1 0]
    @test TDARipserer.adjacency_matrix(LazyWitness(d; nu=3)) == zeros(3, 3)
    @test length(ripserer(LazyWitness(d; nu=3))[1]) == 1
    # A flag octahedral sphere with zero edges has one essential H2 class.
    pairs = [(i, j) for i in 1:6 for j in (i + 1):6 if (i, j) ∉ ((1, 2), (3, 4), (5, 6))]
    octahedral = fill(2.0, length(pairs), 6)
    for (w, (i, j)) in enumerate(pairs)
        octahedral[w, i] = octahedral[w, j] = 0
    end
    h2 = ripserer(LazyWitness(octahedral; nu=0, threshold=1); dim_max=2)[3]
    @test length(h2) == 1
    @test birth(only(h2)) == 0 && !isfinite(only(h2))
    @test_throws ArgumentError LazyWitness(d; nu=4)
    @test_throws ArgumentError Witness(d; max_dimension=-1)

    full = Witness(points, landmarks; max_dimension=3)
    for dimension in 1:3, simplex in full[dimension]
        for facet in TDARipserer.boundary(full, simplex)
            @test birth(facet) <= birth(simplex)
        end
    end
    # Constructor does not silently move omitted high-dimensional killing cells
    # into a truncation intended to retain essential graph cycles.
    @test all(!isfinite, ripserer(Witness(points, landmarks; max_dimension=1))[2])
end
