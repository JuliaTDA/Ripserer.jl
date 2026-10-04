# TDARipserer v0.1.0 (unreleased)

* Add differentiable finite Rips H0/H1 endpoints over F₂, topological losses,
  optional ChainRulesCore/Enzyme extensions, finite-difference tests and executed
  seeded Flux/Lux training examples.
* Add recommended MLJ PH hyperparameter ranges and recursive pipeline tuning.

* Add solver-certified F₂ optimal homologous cycles and filling volumes on
  Rips/Custom through an optional JuMP/HiGHS extension, with deterministic tie handling.
* Add DelayEmbeddings adapters, sliding-window PH and finite-threshold image pipelines.

* Add relaxed weak and lazy witness filtrations, including zero-valued lazy witness edges and landmark provenance.
* Add online sparse F2 zigzag persistence for simplex insertion/deletion streams in arbitrary homology dimensions.
* Add fixed-complex vineyards, incremental boundary-reduction updates, chronological adjacent transpositions and continuous vine identities.

* Establishes an independent experimental JuliaTDA fork of Ripserer.jl with its own package name and UUID.
* Preserves the public function/type names; updates imports, MLJ metadata, tests and documentation to the new package identity.

## Inherited upstream release history

# v0.16.13

Fix performance issues with high-dimensional `Custom` filtrations.

# v0.16.9

Replace LightGraps with Graphs.

# v0.16.4

Add [MLJ.jl](https://github.com/alan-turing-institute/MLJ.jl) support.

# v0.16.3

New function: `midlife`.

# v0.16.0

External changes:

* `progress` keyword argument renamed to `verbose`, `field_type` keyword argument renamed to
  `field`.
* New interface: `ripserer(::Type{AbstractFiltration}, args...; kwargs...)`.
* Added `CircularCoordinates`.

Interface changes:

* All filtration constructors now have to take `verbose` as a keyword argument.
* Replaced vectors of `ChainElement`s with `Chain`s.
* Added `AbstractCell`.
* `Cube` is now an `AbstractCell`, `AbstractSimplex` is reserved for actual simplices.
* Simplices are no longer `Array`s.
* `simplex`, `unsafe_simplex`, and `unsafe_cofacet` no longer take a `sign` argument.

# v0.15.4

* Use `PersistenceDiagrams` v0.8.

# v0.15.3

* Fix type instability in zeroth interval generation.

# v0.15.2

* Update compat with Distances.jl.

# v0.15.1

* Representative cocycles are computed for infinite intervals.

# v0.15.0

* `SparseRips(...)` is deprecated. Use `Rips(...; sparse=true)`.
* `AbstractFiltration`s now need to define `births` instead of `birth`.
* Results are now sorted by persistence instead of birth time.
* Homology is now computed with the `alg=:homology` keyword argument.
* Involuted homology can be computed with the `alg=:involuted` keyword argument.
* The `reps` keyword argument can be set to a collection of integers, finding
  representatives only for specified dimensions.
* Implicit or explicit reduction can be set with the `implicit` keyword argument.
* Improved progress printing.
* New function: `find_apparent_pairs`.
