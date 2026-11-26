# API

## TDARipserer

```@docs
ripserer
```

For exact Euclidean barcodes through degree one, pass coordinates directly to
`ripserer`:

```julia
points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
h0, h1 = ripserer(points; alg=:reduced_rips)
h0, h1 = ripserer(points; alg=:reduced_rips, threshold=1.0)
```

The second call returns the square's loop with an infinite endpoint because the
diagonals lie above the threshold. The backend uses Float64 pointwise distances,
a KD-tree, and the relative neighborhood graph to bound the number of remaining
pairs. The default threshold is the usual Rips radius computed with those
pointwise distances. Finding the relative neighborhood graph at that default
scale can still require a quadratic pair scan. The backend returns birth and
death simplices and supports `dim_max=0` or `1`, modulus 2, and `reps=false`.

For a tight explicit threshold, `ripserer(points; threshold=1.0, radius_graph=true)`
uses the standard cohomology backend with a sparse Euclidean radius graph built
directly from the coordinates. This also supports higher degrees and
representatives. See [the benchmark script](https://github.com/JuliaTDA/Ripserer.jl/blob/master/benchmark/reduced_rips.jl)
for a comparison with the existing Rips paths.

## Filtrations

```@docs
Rips
```

```@docs
Cubical
```

```@docs
Custom
```

```@docs
Alpha
```

```@docs
EdgeCollapsedRips
```

## Persistence Diagrams

Persistence diagrams live in a separate package,
[TDAPersistenceDiagrams.jl](https://github.com/mtsch/PersistenceDiagrams.jl). The package is
documented in detail [here](https://mtsch.github.io/PersistenceDiagrams.jl/dev/).

If you are looking for
Wasserstein or bottleneck distances, persistence images, betti curves, landscapes, and
similar, you will need to run `using TDAPersistenceDiagrams`.

For convenience, the following basic functionality is reexported by TDARipserer:

```@docs
TDAPersistenceDiagrams.PersistenceDiagram
```

```@docs
TDAPersistenceDiagrams.PersistenceInterval
```

```@docs
birth(::TDAPersistenceDiagrams.PersistenceInterval)
```

```@docs
death(::TDAPersistenceDiagrams.PersistenceInterval)
```

```@docs
persistence(::TDAPersistenceDiagrams.PersistenceInterval)
```

```@docs
midlife
```

```@docs
representative(::TDAPersistenceDiagrams.PersistenceInterval)
```

```@docs
birth_simplex(::TDAPersistenceDiagrams.PersistenceInterval)
```

```@docs
death_simplex(::TDAPersistenceDiagrams.PersistenceInterval)
```

```@docs
barcode
```

## Simplices and Representatives

```@docs
TDARipserer.Simplex
```

```@docs
TDARipserer.Cube
```

```@docs
TDARipserer.dim(::TDARipserer.AbstractCell)
```

```@docs
TDARipserer.birth(::TDARipserer.AbstractCell)
```

```@docs
TDARipserer.index(::TDARipserer.AbstractCell)
```

```@docs
TDARipserer.vertices(::TDARipserer.AbstractCell)
```

```@docs
TDARipserer.Chain
```

```@docs
Mod
```

## MLJ.jl Interface

```@docs
TDARipserer.RipsPersistentHomology
```

```@docs
TDARipserer.AlphaPersistentHomology
```

```@docs
TDARipserer.CubicalPersistentHomology
```

## Experimental Features

```@docs
TDARipserer.reconstruct_cycle
```

```@docs
TDARipserer.Partition
```

```@docs
TDARipserer.CircularCoordinates
```

## Abstract Types and Interfaces

```@docs
TDARipserer.AbstractFiltration
```

```@docs
TDARipserer.nv(::TDARipserer.AbstractFiltration)
```

```@docs
TDARipserer.births(::TDARipserer.AbstractFiltration)
```

```@docs
TDARipserer.vertices(::TDARipserer.AbstractFiltration)
```

```@docs
TDARipserer.edges(::TDARipserer.AbstractFiltration)
```

```@docs
TDARipserer.simplex_type
```

```@docs
TDARipserer.simplex
```

```@docs
TDARipserer.unsafe_simplex
```

```@docs
TDARipserer.unsafe_cofacet
```

```@docs
TDARipserer.threshold(::TDARipserer.AbstractFiltration)
```

```@docs
TDARipserer.columns_to_reduce
```

```@docs
TDARipserer.emergent_pairs
```

```@docs
TDARipserer.postprocess_diagram
```

```@docs
TDARipserer.distance_matrix
```

```@docs
TDARipserer.AbstractRipsFiltration
```

```@docs
TDARipserer.adjacency_matrix(::TDARipserer.AbstractFiltration)
```

```@docs
TDARipserer.AbstractCustomFiltration
```

```@docs
TDARipserer.simplex_dicts
```

```@docs
TDARipserer.AbstractCell
```

```@docs
TDARipserer.AbstractSimplex
```

```@docs
Base.sign(::TDARipserer.AbstractCell)
```

```@docs
Base.:-(::TDARipserer.AbstractCell)
```

```@docs
TDARipserer.coboundary
```

```@docs
TDARipserer.boundary
```
