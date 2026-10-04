# Witness complexes and dynamic persistence

## Witness filtrations

Select landmarks from a large point cloud, then use all points as witnesses.
`Witness` implements a relaxed weak witness complex: every face must have a
weak witness. `LazyWitness` uses witnessed edge values and their flag closure,
which permits implicit PH without materializing every higher simplex. These
are different constructions; neither promises the exact original Rips barcode.

Matrix inputs have witnesses in rows and landmarks in columns. Filtration
vertex labels are positions in the landmark list, and `flt.landmarks` maps them
back to original data indices. The exact weak witness constructor has a
combinatorial landmark cost; limit `max_dimension` and use a threshold.

```julia
points = [[cos(t), sin(t)] for t in range(0, 2π; length=513)[1:512]]
landmarks = collect(1:32:512)
diagrams = ripserer(LazyWitness(points, landmarks; threshold=0.8))
```

## Zigzag persistence

Insert and remove simplices using an online `ZigzagPersistence` state. This
maintains actual homological maps under inclusions, rather than matching
separate diagrams. Coefficients are F2; arbitrary requested homology dimensions
are supported. A deletion must remove a maximal simplex, and an insertion must
follow its faces. Invalid events are rejected before advancing the stream.

State index zero is the initial complex. Each event advances the index by one.
An interval `[b,d)` is active at indices `b` through `d-1`. Infinite deaths are
right-censored by the stream end. Initial maximal simplices are expanded to
their face closure at index zero.

```julia
state = ZigzagPersistence(initial=[(1, 2, 3)])
push!(state, :remove => (1, 2, 3))
push!(state, :remove => (1, 2))
zigzag_diagrams(state)[2] # H1 is [1,2)
```

The sparse bases retain cycles, boundaries, and filling chains; operations
only change the affected chain degrees. The right-filtration ordering
determines which generator closes. The independent test oracle constructs
homology from dense boundary matrices and compares the limit-to-colimit rank
of every subinterval. This distinguishes correct interval pairings from merely
correct Betti numbers.

## Vineyards and incremental persistence

`Vineyard` retains the factorization `R = D*V` on a fixed simplicial complex.
`update_vineyard!` follows affine changes of its filtration values. Each order
crossing conjugates and repairs the existing factorization through an adjacent
transposition. Order-preserving changes require no column additions.

Vine identity follows the unchanged endpoint when birth or death coordinates
cross; it is not obtained by matching newly recomputed diagrams. The knot
sequences in `state.vines` include diagonal crossings. At simultaneous ties,
continuation can be nonunique; the implementation chooses the leftmost
available swap deterministically. Updates validate the final face inequalities
before mutating the state; affine interpolation then also satisfies them.

Supply partial simplex-value updates or a full value vector in the current
`state.simplices` order. A large permutation can involve quadratically many
crossings and may be slower than computing one final diagram. Use zigzag for
changes to the complex itself.

Run `examples/p3_algorithms.jl` for all three workflows, and `benchmark/p3.jl`
for reproducible timings, allocation measurements and operation counters.

References:

- de Silva and Carlsson, [Topological estimation using witness complexes](https://doi.org/10.2312/SPBG/SPBG04/157-166).
- Carlsson and de Silva, [Zigzag Persistence](https://arxiv.org/abs/0812.0197).
- Cohen-Steiner, Edelsbrunner and Morozov, [Vines and vineyards](https://mrzv.org/publications/vineyards/).

```@docs
Witness
LazyWitness
ZigzagEvent
ZigzagPersistence
zigzag_persistence
zigzag_diagrams
VinePoint
Vineyard
update_vineyard!
vineyard_diagrams
```
