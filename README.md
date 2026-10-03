# TDARipserer.jl

Experimental persistent homology for the [JuliaTDA](https://github.com/JuliaTDA)
ecosystem, derived from [Ripserer.jl](https://github.com/mtsch/Ripserer.jl)
by Matija Čufar and contributors. The fork is maintained by G. Vituri and
JuliaTDA contributors and starts its own version series at `0.1.0`.

Package UUID: `d8d97730-eb17-4b7a-9fb6-200f7d8ece7d`.
The GitHub repository currently retains its original `JuliaTDA/Ripserer.jl`
slug; the Julia package and local checkout are named `TDARipserer` and
`TDARipserer.jl`. This fork is not registered.

[![Build Status](https://github.com/JuliaTDA/Ripserer.jl/actions/workflows/Test.yml/badge.svg?branch=master)](https://github.com/JuliaTDA/Ripserer.jl/actions/workflows/Test.yml)

## Fork scope and compatibility

This is an independent home for experimental filtrations, algorithms and
JuliaTDA integration. Its starting point is upstream Ripserer `0.16.14`
(commit `9583a3b`). No new homology algorithm is introduced by the rename.
The functions and types such as `ripserer`, `Rips` and `Cubical` keep their names.

TDARipserer depends on **TDAPersistenceDiagrams**, so its output uses that
fork's `PersistenceDiagram` type. Types from the original packages and these
forks are distinct; downstream code must import the corresponding fork.
API compatibility is experimental and breaking changes increment the minor
version before `1.0`. See [CHANGELOG.md](CHANGELOG.md).

## Tracking upstream

The Git history and upstream citations are preserved. Inherited Git tags
describe upstream releases, rather than releases of this independent fork. A local `upstream`
remote points to `https://github.com/mtsch/Ripserer.jl.git`. Fresh clones can
add it with `git remote add upstream https://github.com/mtsch/Ripserer.jl.git`.
Use `git fetch --no-tags upstream` to inspect changes, then selectively merge
or cherry-pick relevant fixes. Preserve the fork's name, UUID, dependency
identity and MLJ metadata, and port changes to the renamed entry point
`src/TDARipserer.jl`. Run tests against TDAPersistenceDiagrams after syncing.

## Introduction

TDARipserer is a pure Julia implementation of the [Ripser](https://github.com/Ripser/ripser)
algorithm for computing persistent homology. Its aims are to be easy to use, generic, and
fast.

See the [documentation sources](docs/src/index.md) for more information and
usage examples.

If you're looking for persistence diagram-related functionality such as Wasserstein or
bottleneck distances, persistence images, or persistence curves, please see
[TDAPersistenceDiagrams.jl](https://github.com/JuliaTDA/PersistenceDiagrams.jl).

## Quick start

For local development, activate your working environment and develop the
two sibling checkouts together (paths below are relative to that environment):

```julia
using Pkg
Pkg.develop([
    PackageSpec(path = "../TDAPersistenceDiagrams.jl"),
    PackageSpec(path = "../TDARipserer.jl"),
])
using TDARipserer
```

For the repository URLs, use a batch `Pkg.develop` with
`https://github.com/JuliaTDA/PersistenceDiagrams.jl` and
`https://github.com/JuliaTDA/Ripserer.jl` after these changes are pushed.

Now, generate some data.

```julia
julia> data = [(rand(), rand(), rand()) for _ in 1:200]
```

The main exported function in this package is
[`ripserer`](https://mtsch.github.io/Ripserer.jl/dev/api/ripserer/#Ripserer.ripserer). By
default, it computes Vietoris-Rips persistent homology on point cloud data and distance
matrices.

```julia
julia> ripserer(data)
# 2-element Vector{TDAPersistenceDiagrams.PersistenceDiagram}:
#  200-element 0-dimensional PersistenceDiagram
#  84-element 1-dimensional PersistenceDiagram
```

[Several other filtration
types](https://mtsch.github.io/Ripserer.jl/dev/api/ripserer/#Filtrations) are supported. We
tell `ripserer` to use them by passing them as the first argument.

```julia
julia> ripserer(EdgeCollapsedRips, data)
# 2-element Vector{TDAPersistenceDiagrams.PersistenceDiagram}:
#  200-element 0-dimensional PersistenceDiagram
#  84-element 1-dimensional PersistenceDiagram
```

Sometimes you may want to initialize a filtration in advance.

```julia
julia> rips = EdgeCollapsedRips(data, threshold=1)
# EdgeCollapsedRips{Int64, Float64}(nv=200)
```
```julia
julia> ripserer(rips, dim_max=2)
# 3-element Vector{TDAPersistenceDiagrams.PersistenceDiagram}:
#  200-element 0-dimensional PersistenceDiagram
#  84-element 1-dimensional PersistenceDiagram
#  16-element 2-dimensional PersistenceDiagram
```

TDARipserer supports plotting with
[Plots.jl](https://github.com/JuliaPlots/Plots.jl). Makie plotting for the fork
is provided by [TDAplots.jl](https://github.com/JuliaTDA/TDAplots.jl).

Plotting persistence diagrams and barcodes is straightforward:

```julia
using Plots
result = ripserer(data, dim_max=2)
plot(plot(result), barcode(result))
```
![](docs/src/assets/readme-plot-1.svg)

```julia
barcode(result)
```
![](docs/src/assets/readme-plot-2.svg)

## Attribution and license

The original MIT copyright and license are retained in [LICENSE](LICENSE),
with an additional copyright notice for JuliaTDA contributions. For the
underlying Ripserer algorithms, retain the upstream citations in
[CITATION.bib](CITATION.bib); those publications describe upstream work.
