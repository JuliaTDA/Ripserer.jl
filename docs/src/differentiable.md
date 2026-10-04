# Differentiable Rips persistence

`differentiable_persistence(X)` returns the finite diagram as a matrix of birth
and death coordinates. Columns of `X` are Euclidean input points. Alternatively,
`input=:distances` accepts a dense symmetric nonnegative distance matrix with
zero diagonal, giving gradients with respect to filtration edge values.

The supported subset is **Rips H0/H1 over F₂**, with the ordinary dense reduction,
without cutoff, threshold, collapse, or sparsification. Essential bars are
excluded. Gradients follow the critical edges of birth/death simplices. They
are exact on a fixed pairing stratum; ties and pairing changes are nonsmooth.
At ties, the lexicographically first maximizing edge selects a branch
subgradient. Coincident points use the zero subgradient. Recompute the pairing
after every optimizer step. No derivative through the discrete pairing itself
is claimed. Distance gradients split equally across symmetric entries.

```julia
using TDARipserer, Zygote
X = [0.0 1.0 1.23 0.11; 0.0 0.07 1.1 1.02]
diagram = differentiable_persistence(X;dim=1)
g = Zygote.gradient(x -> topological_loss(x;target=0.7),X)[1]
```

`topological_loss` sums finite lifetimes raised to `power` (default 2).
`target=...` takes the squared discrepancy from a target total;
`negate=true` changes the sign for maximization. Maximizing persistence alone
can expand scale: use an application-specific data or regularization loss.
Empty diagrams have zero loss and zero coordinate gradient.

Julia ≥1.9 loads optional AD extensions when `Zygote`/`ChainRulesCore` or
`Enzyme` is loaded. With Enzyme:

```julia
using Enzyme
g = zeros(size(X))
Enzyme.autodiff(Enzyme.Reverse, x -> topological_loss(x;target=0.7),
    Enzyme.Active, Enzyme.Duplicated(X,g))
```

Only pairing computation is marked inactive in Enzyme; numeric endpoint
extraction remains active. The two-argument form with `PersistencePairing`
holds a supplied pairing fixed and is useful for local calculations.

## Flux and Lux training

Runnable seeded examples are in `examples/differentiable/flux.jl` and
`examples/differentiable/lux.jl`. Both train a small linear map on a noisy circle
to match a target H1 total persistence and assert that the objective decreases.
They use Julia ≥1.10, Flux 0.16, Lux 1, Zygote 0.7, and Optimisers 0.4.
From the ecosystem workspace:

```sh
julia --project=TDARipserer.jl/examples/differentiable -e 'using Pkg; Pkg.develop([PackageSpec(path="TDAPersistenceDiagrams.jl"),PackageSpec(path="TDARipserer.jl")]); Pkg.instantiate()'
julia --project=TDARipserer.jl/examples/differentiable TDARipserer.jl/examples/differentiable/flux.jl
julia --project=TDARipserer.jl/examples/differentiable TDARipserer.jl/examples/differentiable/lux.jl
julia --project=TDARipserer.jl/examples/differentiable TDARipserer.jl/test/extra/differentiable.jl
```

Tests compare births and deaths separately against central finite differences,
check translation invariance and symmetric-distance directions, and verify
Zygote and Enzyme reverse modes. AD and both training examples have a dedicated
CI job. These examples demonstrate the supported layer, not general convergence.

Reference: [Carrière et al. (2021), Optimizing persistent homology based
functions](https://proceedings.mlr.press/v139/carriere21a.html).

```@docs
PersistencePairing
persistence_pairing
differentiable_persistence
topological_loss
```
