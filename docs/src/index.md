# TDARipserer.jl

_Flexible and efficient persistent homology computation._

An independent experimental JuliaTDA fork of
[Ripserer.jl](https://github.com/mtsch/Ripserer.jl), originally authored by
Matija Čufar ([@mtsch](https://github.com/mtsch/)) and contributors.
Fork maintained by G. Vituri and JuliaTDA contributors.
The fork has its own UUID and version series; see its [README](https://github.com/JuliaTDA/Ripserer.jl).

## Introduction

TDARipserer is a pure Julia library for computing persistent homology based on the
[Ripser](https://github.com/Ripser/ripser) algorithm. Roughly speaking, persistent homology
detects the global topological and local geometric structure of data in a noise-resistant,
stable way. If you are unfamiliar with persistent homology, I recommend reading this
[excellent
introduction](https://towardsdatascience.com/persistent-homology-with-examples-1974d4b9c3d0).

Please see the [Usage Guide](@ref) for a quick introduction, and the [API](@ref) page for
detailed descriptions of TDARipserer's functionality.

While this package is fully functional, it is still in development and should not be
considered stable. I try to disrupt the public interface as little as possible, but breaking
changes might still occur from time to time.

## Installation

This fork is not registered. For local development from a sibling package
environment, develop both checkouts together:

```julia
using Pkg
Pkg.develop([
    PackageSpec(path = "../TDAPersistenceDiagrams.jl"),
    PackageSpec(path = "../TDARipserer.jl"),
])
using TDARipserer
```

Julia 1.6 or later is supported. The public functions keep their upstream
names, but types from the forks are distinct from the original packages.

## Features

TDARipserer and its companion package
[TDAPersistenceDiagrams.jl](https://github.com/JuliaTDA/PersistenceDiagrams.jl) currently support

* Fast Vietoris-Rips and cubical, and alpha complex persistent homology computation.
* Representative cocycle, cycle, and critical simplex computation.
* Convenient persistence diagram and representative cocycle visualization via
  [Plots.jl](https://github.com/JuliaPlots/Plots.jl). Makie plotting for the fork
  is provided by [TDAplots.jl](https://github.com/JuliaTDA/TDAplots.jl).
* Bottleneck and Wasserstein matching and distance computation.
* Various persistence diagram vectorization functions, implemented with persistence images
  and persistence curves.
* Easy extensibility through a documented API.
* Integration with [MLJ.jl](https://github.com/alan-turing-institute/MLJ.jl).
* Experimental shortest representative cycle computation.
* Experimental sparse circular coordinate computation.

To access some of the features, you need to use the TDAPersistenceDiagrams.jl package.

## Performance

Much like Ripser, TDARipserer uses several computational tricks to achieve its speed. Among
others, these include an implicit simplicial complex representation and the clearing
optimization. For a more detailed overview of these optimizations, check out [Ulrich Bauer's
article on Ripser](https://arxiv.org/abs/1908.02518).

In general, the performance of TDARipserer is very close to
[Ripser](https://github.com/Ripser/ripser), usually within around 30%. TDARipserer's strength
performance-wise is very sparse inputs, where it can sometimes outperform Ripser. It also
computes some things Ripser skips, like the critical simplices.

TDARipserer's Cubical homology is up to 3× slower than that of [Cubical
Ripser](https://github.com/CubicalRipser/), which uses a more specialized
algorithm. TDARipserer is still a good choice for small 3d images and large 2d images. Unlike
Cubical Ripser, it also supports computations on images of dimensions higher than 4.

See the [Benchmarks](@ref) section for more detailed benchmarks.

## Extending

TDARipserer is designed to be easily extended with new simplex or filtration types. See the
[Abstract Types and Interfaces](@ref) API section for more information.

If you have written an extension or are having trouble implementing one, please feel free to
open a pull request or an issue in the [JuliaTDA fork](https://github.com/JuliaTDA/Ripserer.jl).

## Contributing

Contributions and experimental ideas are welcome in the
[JuliaTDA fork](https://github.com/JuliaTDA/Ripserer.jl).
Upstream improvements can be ported selectively while preserving this fork's
name, UUID and dependency identities. The original MIT license and credits
are retained.

## Citing

The underlying upstream Ripserer implementation is described in the [JOSS
paper](https://joss.theoj.org/papers/10.21105/joss.02614).

A bibtex entry is provided in
[CITATION.bib](https://github.com/mtsch/Ripserer.jl/blob/master/CITATION.bib).
