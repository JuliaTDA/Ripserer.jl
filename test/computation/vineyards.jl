using TDARipserer
using TDAPersistenceDiagrams
using Test
using Random

endpoints(diagrams) = [sort([(birth(x), death(x)) for x in d]) for d in diagrams]

function verify_factorization(state)
    n = length(state.simplices)
    positions = Dict(s => i for (i, s) in enumerate(state.simplices))
    d = [BitSet([positions[f] for f in TDARipserer._facets(s)]) for s in state.simplices]
    for j in 1:n
        result = BitSet()
        for i in state.basis[j]
            symdiff!(result, d[i])
        end
        @test result == state.reduced[j]
        @test all(i -> i <= j, state.basis[j])
        @test j in state.basis[j]
    end
    lows = [maximum(c) for c in state.reduced if !isempty(c)]
    @test allunique(lows)
end

@testset "Incremental persistence and continuous vines" begin
    simplices = TDARipserer._simplex_closure([(1, 2, 3, 4)])
    rng = MersenneTwister(31)
    initial_values = Dict(s => 2 * (length(s) - 1) + rand(rng) for s in simplices)
    state = Vineyard(collect(initial_values); dim_max=2)
    ids = Set(keys(state.vines))
    for frame in 1:12
        target = Dict(s => 2 * (length(s) - 1) + rand(rng) for s in simplices)
        update_vineyard!(state, target; time=frame)
        expected = ripserer(Custom(collect(target)); dim_max=2)
        @test endpoints(vineyard_diagrams(state)) == endpoints(expected)
        @test Set(keys(state.vines)) == ids
        verify_factorization(state)
    end
    @test state.transpositions > 0
    @test issorted([swap.time for swap in state.swaps])
    # At every crossing, a vine may change its defining simplices but cannot
    # jump in diagram coordinates. Consecutive knots at equal times coincide.
    for knots in values(state.vines), i in 2:length(knots)
        if knots[i].time == knots[i - 1].time
            @test knots[i].birth ≈ knots[i - 1].birth
            @test knots[i].death ≈ knots[i - 1].death
        end
    end
    additions = state.column_additions
    update_vineyard!(state, Dict(s => v + 1 for (s, v) in state.values))
    @test state.column_additions == additions
    @test_throws ArgumentError update_vineyard!(state, Dict((1, 2) => -1))
    @test_throws ArgumentError update_vineyard!(state, Dict((5,) => 0))
    @test_throws ArgumentError update_vineyard!(state, Dict(); time=0)
    @test_throws ArgumentError Vineyard([(1,) => 2, (1, 2) => 1])
    @test endpoints(vineyard_diagrams(Vineyard(Custom([(1, 2, 3) => 2, (1, 2) => 1])))) ==
        endpoints(ripserer(Custom([(1, 2, 3) => 2, (1, 2) => 1])))

    # Cell dimensions need not occupy separate value ranges. These lower-star
    # frames contain simultaneous ties and crossings between unrelated cells
    # of different dimensions, including a completely collapsed filtration.
    tied = Vineyard([s => 0.0 for s in simplices]; dim_max=2)
    for frame in 1:6
        vertex_values = 8 .* rand(rng, 4)
        target = Dict(s => maximum(vertex_values[collect(s)]) for s in simplices)
        update_vineyard!(tied, target; time=frame)
        @test endpoints(vineyard_diagrams(tied)) ==
            endpoints(ripserer(Custom(collect(target)); dim_max=2))
        verify_factorization(tied)
    end
    update_vineyard!(tied, zeros(length(tied.simplices)); time=7)
    @test endpoints(vineyard_diagrams(tied)) ==
        endpoints(ripserer(Custom([s => 0.0 for s in simplices]); dim_max=2))
    verify_factorization(tied)
    for knots in values(tied.vines), i in 2:length(knots)
        if knots[i].time == knots[i - 1].time
            @test knots[i].birth ≈ knots[i - 1].birth
            @test knots[i].death ≈ knots[i - 1].death
        end
    end
end
