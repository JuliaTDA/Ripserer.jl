# Optimal chains, time series and images

## Explanatory cycles and filling volumes

On Julia 1.9 or newer, load JuMP and HiGHS, then compute an interval with representatives:

```julia
using TDARipserer, JuMP, HiGHS
f = Rips([[0.,0.], [1.,0.], [1.,1.], [0.,1.]])
interval = only(ripserer(f; reps=true)[2])
cycle = optimal_cycle(f, interval)       # minimum number of edges
volume = optimal_volume(f, interval)    # minimum number of filling triangles
```

The optimization is a binary/integer F₂ program, not the oriented integer/real LP
of [Dey, Hirani and Krishnamoorthy](https://arxiv.org/abs/1001.0338). A homologous chain
satisfies z=c+∂y mod 2; a filling satisfies ∂z=c mod 2. Nonnegative simplex weights
can be passed as a function. Returned `filling` certifies the cycle's homology;
optimization statuses and parity certificates are checked independently.

Rips and Custom simplicial filtrations are supported. H₁ cocycles are converted
with `reconstruct_cycle`. For `alg=:homology` representatives, including higher
dimensions, pass `representative_kind=:cycle`. Other coefficient fields, cubical
cells and arbitrary cocycle-to-cycle conversion are not supported.

Cycle threshold `at` must lie in the selected interval; the interval volume helper
reconstructs at birth and fills at death. Infinite deaths require an explicit finite
filling threshold. Nonbounding cycles raise an infeasibility error. `max_cells`
limits enumeration, which can be expensive on dense complexes. Optional `time_limit`
limits each solver invocation; a timeout does not return an unproven optimum.

Default tie handling minimizes the support vector lexicographically among chains
with cost at most optimum plus `tolerance=1e-8`. This takes one extra solve per cell.
Use `tie_break=:solver` to obtain any optimal support efficiently. No uniqueness of
an explanatory geometric chain is implied.

## Time series and cubical inputs

On Julia 1.9+, `delay_embedding` delegates to DelayEmbeddings.embed; `suggest_delay` delegates to
its estimator. `sliding_window_ph` consumes already embedded points, including native
StateSpaceSets from DynamicalSystems, and records each embedded window's indices.
Delay selection is exploratory; use chronological or blocked validation. Embeddings
lose `(dimension-1)*delay` samples, and overlapping windows are dependent. Two runnable
examples are in JuliaTDA.jl/examples/p2, covering a periodic signal and a Hénon trajectory.

`image_persistence(image; transform=identity)` accepts callable preprocessing, such
as MetricSpaces.ImageFiltration. It chooses the maximum finite value as its default
Cubical threshold, excluding Inf background vertices. An all-excluded image is an
error. The convention is vertex lower-star, matching MetricSpaces.cubical_ecc.

```@docs
optimal_cycle
optimal_volume
OptimalChainResult
delay_embedding
suggest_delay
sliding_window_ph
image_persistence
```
