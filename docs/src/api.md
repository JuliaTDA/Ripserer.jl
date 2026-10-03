# API

## TDARipserer

```@docs
ripserer
```

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
