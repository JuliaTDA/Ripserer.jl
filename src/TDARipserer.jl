"""
# TDARipserer.jl

Efficient computation of persistent homology.

Experimental JuliaTDA fork of Ripserer.jl. See the documentation in `docs/src`.
"""
module TDARipserer

using Compat

using Base: @propagate_inbounds
using Base.Cartesian
using LinearAlgebra
using SparseArrays

using DataStructures
using Distances
using IterTools
using Graphs
using MiniQhull
using TDAPersistenceDiagrams
using ProgressMeter
using RecipesBase
using StaticArrays
using TupleTools

import MLJModelInterface

# This functionality is imported to avoid having to deal with name clashes. There is no
# piracy involved here.
import Graphs: vertices, edges, nv, adjacency_matrix

# Reexporting basics makes TDARipserer usable without having to import another package.
import TDAPersistenceDiagrams: birth, threshold, dim
export birth,
    death, persistence, midlife, representative, birth_simplex, death_simplex, barcode

export Mod,
    Simplex,
    Cube,
    index,
    vertices,
    dim,
    threshold,
    simplex,
    coefficient,
    Rips,
    Cubical,
    Custom,
    Alpha,
    EdgeCollapsedRips,
    Witness,
    LazyWitness,
    ZigzagEvent,
    ZigzagPersistence,
    zigzag_persistence,
    zigzag_diagrams,
    VinePoint,
    Vineyard,
    update_vineyard!,
    vineyard_diagrams,
    ripserer,
    reconstruct_cycle,
    Partition,
    CircularCoordinates,
    RipsPersistentHomology,
    AlphaPersistentHomology,
    CubicalPersistentHomology

include("base/primefield.jl")
include("base/abstractcell.jl")
include("base/abstractfiltration.jl")
include("base/abstractsimplex.jl")
include("base/simplex.jl")
include("base/chainelement.jl")
include("base/chain.jl")
include("base/simplexrecipes.jl")

include("computation/utils.jl")
include("computation/zerodimensional.jl")
include("computation/reducedmatrix.jl")
include("computation/cohomology.jl")
include("computation/ripserer.jl")

include("filtrations/utils.jl")
include("filtrations/rips.jl")
include("filtrations/cubical.jl")
include("filtrations/custom.jl")
include("filtrations/alpha.jl")
include("filtrations/edgecollapse.jl")
include("filtrations/witness.jl")
include("computation/binary_algebra.jl")
include("computation/zigzag.jl")
include("computation/vineyards.jl")

include("extra/cycles.jl")
include("extra/circularcoordinates.jl")
include("extra/mlj.jl")
include("extra/differentiable.jl")
export PersistencePairing, persistence_pairing, differentiable_persistence, topological_loss

include("extra/optimal_chains.jl")
export OptimalChainResult, optimal_cycle, optimal_volume

include("extra/workflows.jl")
export delay_embedding, suggest_delay, sliding_window_ph, image_persistence

end
