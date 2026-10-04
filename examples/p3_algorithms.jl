using TDARipserer
using TDAPersistenceDiagrams

# P3.1: 512 witnesses, 16 landmarks, implicit PH on a lazy witness complex.
points = [[cos(t), sin(t)] for t in range(0, 2π; length=513)[1:512]]
landmarks = collect(1:32:512)
lazy = LazyWitness(points, landmarks; nu=1, threshold=0.8)
witness_diagrams = ripserer(lazy)
@assert !isempty(witness_diagrams[2])
println("Lazy witness H1: ", [(birth(x), death(x)) for x in witness_diagrams[2]])
full = Witness(points, landmarks; max_dimension=2, threshold=0.8)
@assert !isempty(ripserer(full)[2])

# P3.2: a filled triangle loses its face, gains a loop, then loses an edge.
stream = ZigzagPersistence(; initial=[(1, 2, 3)])
push!(stream, :remove => (1, 2, 3))
push!(stream, :remove => (1, 2))
@assert [(birth(x), death(x)) for x in zigzag_diagrams(stream)[2]] == [(1, 2)]
println("Zigzag H1: ", zigzag_diagrams(stream)[2])

# P3.4: the same triangle with a moving filtration, retaining the reduction.
filtration = [(1,) => 0.0, (2,) => 0.0, (3,) => 0.0,
              (1, 2) => 1.0, (2, 3) => 2.0, (1, 3) => 3.0, (1, 2, 3) => 4.0]
state = Vineyard(filtration)
update_vineyard!(state, Dict((1, 2) => 3.0, (1, 3) => 1.0); time=1.0)
println("Vineyard H1: ", vineyard_diagrams(state)[2])
println("Adjacent swaps: ", state.transpositions, "; column additions: ", state.column_additions)
@assert state.transpositions > 0
