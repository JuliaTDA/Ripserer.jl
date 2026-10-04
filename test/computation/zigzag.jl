using TDARipserer
using TDAPersistenceDiagrams
using Test
using Random
include("zigzag_oracle.jl")

@testset "Online zigzag persistence" begin
    events = [:add => (1,), :add => (2,), :add => (1, 2), :remove => (1, 2)]
    d = only(zigzag_persistence(events; dim_max=0))
    @test sort([(birth(x), death(x)) for x in d]) == [(1, Inf), (2, 3), (4, Inf)]
    state = ZigzagPersistence(; initial=[(1, 2, 3)], dim_max=2)
    push!(state, :remove => (1, 2, 3))
    @test birth(only(zigzag_diagrams(state)[2])) == 1
    push!(state, :remove => (1, 2))
    @test (birth(only(zigzag_diagrams(state)[2])), death(only(zigzag_diagrams(state)[2]))) == (1, 2)
    @test_throws ArgumentError push!(state, :remove => (1,))
    @test_throws ArgumentError push!(state, :add => (1,))
    @test_throws ArgumentError push!(state, :add => (1, 2, 3))
    @test state.steps == 2
    @test_throws ArgumentError ZigzagEvent(:wrong, (1,))
    @test_throws ArgumentError ZigzagEvent(:add, (1, 1))
    stream = ZigzagPersistence(initial=[(1, 2), (2, 3), (1, 3)])
    for i in 1:2000
        push!(stream, ZigzagEvent(isodd(i) ? :add : :remove, (1, 2, 3)))
    end
    @test stream.next_id - 1 == stream.max_active_simplices == 7
    @test length(zigzag_diagrams(stream)[2]) == 1001

    # S^2 -> filled tetrahedron -> S^2: dimensions above H1 are supported.
    sphere = [(1, 2, 3), (1, 2, 4), (1, 3, 4), (2, 3, 4)]
    state = ZigzagPersistence(; initial=sphere, dim_max=2)
    @test birth(only(zigzag_diagrams(state)[3])) == 0
    push!(state, :add => (1, 2, 3, 4))
    push!(state, :remove => (1, 2, 3, 4))
    @test sort([(birth(x), death(x)) for x in zigzag_diagrams(state)[3]]) == [(0, 1), (2, Inf)]

    # Compare all subinterval ranks to a dense limit-to-colimit construction
    # on independently computed homology spaces, not just Betti numbers.
    for seed in 1:4
        rng = MersenneTwister(seed)
        universe = TDARipserer._simplex_closure([(1, 2, 3, 4)])
        active = Set{Tuple}(seed <= 2 ? universe[1:(end - 1)] : Tuple[])
        events = ZigzagEvent[]
        snapshots = [copy(active)]
        for _ in 1:22
            additions = [s for s in universe if !(s in active) &&
                all(f -> f in active, TDARipserer._facets(s))]
            removals = [s for s in active if !any(t -> length(t) > length(s) &&
                all(v -> v in t, s), active)]
            choices = vcat([ZigzagEvent(:add, s) for s in additions],
                           [ZigzagEvent(:remove, s) for s in sort(removals)])
            event = rand(rng, choices)
            push!(events, event)
            event.operation == :add ? push!(active, event.simplex) : delete!(active, event.simplex)
            push!(snapshots, copy(active))
        end
        diagrams = zigzag_persistence(events; dim_max=2, initial=collect(first(snapshots)))
        directions = [e.operation == :add for e in events]
        for dimension in 0:2
            spaces = [oracle_homology(c, dimension) for c in snapshots]
            maps = [directions[i] ? oracle_inclusion(spaces[i], spaces[i + 1]) :
                oracle_inclusion(spaces[i + 1], spaces[i]) for i in eachindex(events)]
            for a in eachindex(snapshots), b in a:length(snapshots)
                expected = oracle_generalized_rank(spaces, maps, directions, a, b)
                actual = count(x -> birth(x) <= a - 1 && death(x) > b - 1, diagrams[dimension + 1])
                @test actual == expected
            end
        end
    end
end
