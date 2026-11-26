struct LazyFalses end
Base.getindex(::LazyFalses, i, j) = false

mutable struct _EdgeMarks
    positions::Dict{Tuple{Int,Int},Int}
    marks::Vector{Int}
    generation::Int
end

_EdgeMarks(positions) = _EdgeMarks(positions, zeros(Int, length(positions)), 1)
_reset!(marks::_EdgeMarks) = (marks.generation += 1)

function Base.getindex(marks::_EdgeMarks, u, v)
    i = get(marks.positions, minmax(u, v), 0)
    return i != 0 && marks.marks[i] == marks.generation
end

function _mark!(marks::_EdgeMarks, u, v)
    u == v && return
    marks.marks[marks.positions[minmax(u, v)]] = marks.generation
    return nothing
end

function _get_edge_neighbors!(buffer, graph::AbstractGraph, edge, removed=LazyFalses())
    s, d = src(edge), dst(edge)
    nu = neighbors(graph, s)
    nv = neighbors(graph, d)
    empty!(buffer)
    # Sorted version of intersect
    i, j = 1, 1
    while i ≤ length(nu) && j ≤ length(nv)
        u = nu[i]
        v = nv[j]
        if u == v
            if !removed[u, s] && !removed[u, d]
                push!(buffer, u)
            end
            i += 1
            j += 1
        elseif u < v
            i += 1
        elseif v < u
            j += 1
        end
    end
    return buffer
end

function _get_adjacent_edges!(adjacent, graph::AbstractGraph, edge, removed=LazyFalses())
    for u in (src(edge), dst(edge))
        for v in neighbors(graph, u)
            if !removed[v, u]
                _mark!(adjacent, u, v)
            end
        end
    end
    return adjacent
end

function _not_dominated(buffer, graph::AbstractGraph, edge, removed=LazyFalses())
    # N[e] ⊆ N[v] for some v ∉ edge
    s, d = src(edge), dst(edge)
    ne = _get_edge_neighbors!(buffer, graph, edge, removed)
    for w in ne
        (w == s || w == d) && continue

        nw = neighbors(graph, w)
        # Sorted version of issubset
        i, j = 1, 1
        while i ≤ length(ne) && j ≤ length(nw)
            u = ne[i]
            v = nw[j]
            if removed[v, w]
                #
            elseif u < v
                break
            elseif u == v
                i += 1
            end
            j += 1
        end
        i > length(ne) && return false
    end
    return true
end

function _core_graph(filtration::Rips{<:Any,T}; verbose=false) where {T}
    simplices = sort!(edges(filtration))
    endpoints = [minmax(Int.(vertices(e))...) for e in simplices]
    positions = Dict(pair => i for (i, pair) in enumerate(endpoints))

    buffer = Int[]
    removed = _EdgeMarks(positions)
    adjacent = _EdgeMarks(positions)
    retained = falses(length(simplices))
    weights = Vector{T}(undef, length(simplices))
    graph = Graph(nv(filtration))
    # Add self loops so i ∈ neighbors(graph, i)
    for i in 1:nv(filtration)
        add_edge!(graph, i, i)
    end

    if verbose
        progbar = Progress(length(simplices); desc="Collapsing edges...")
    end
    for (i, (u, v)) in enumerate(endpoints)
        curr_edge = Edge(u, v)
        curr_time = birth(simplices[i])
        add_edge!(graph, curr_edge)
        if _not_dominated(buffer, graph, curr_edge)
            retained[i] = true
            weights[i] = curr_time

            _reset!(adjacent)
            _reset!(removed)
            _get_adjacent_edges!(adjacent, graph, curr_edge)
            # Some edges that were not dominated before may be dominated now.
            for j in (i - 1):-1:1
                retained[j] && continue
                edge = Edge(endpoints[j]...)
                if adjacent.marks[j] == adjacent.generation &&
                        _not_dominated(buffer, graph, edge, removed)
                    retained[j] = true
                    weights[j] = curr_time
                    _get_adjacent_edges!(adjacent, graph, edge, removed)
                else
                    removed.marks[j] = removed.generation
                end
            end
        end
        verbose && next!(progbar)
    end
    rows, cols, values = Int[], Int[], T[]
    for i in eachindex(retained)
        retained[i] || continue
        u, v = endpoints[i]
        push!(rows, u, v)
        push!(cols, v, u)
        push!(values, weights[i], weights[i])
    end
    for (v, value) in enumerate(births(filtration))
        iszero(value) && continue
        push!(rows, v)
        push!(cols, v)
        push!(values, value)
    end
    return sparse(rows, cols, values, nv(filtration), nv(filtration))
end

"""
    EdgeCollapsedRips{I, T} <: AbstractRipsFiltration{I, T}

Perform a sequence of edge collapses on a filtration. This may significantly reduce
computation time and does not change the result. The speedup is especially apparent with
datasets that have a boundary, and with high-dimensional persistent homology computation.

Doing the collapses can be time-consuming. The construction uses storage linear in the
number of vertices and edges of the input filtration.

See the reference below for a description of the algorithm.

# Constructors

* `EdgeCollapsedRips(::Rips; verbose=false, threshold=nothing)`:
  Collapse a given filtration, retaining its threshold by default. Setting `verbose` shows
  a progress bar.

* `EdgeCollapsedRips(::EdgeCollapsedRips; verbose=false, threshold=nothing)`: Allows
  changing `I` or `threshold` without recomputing.

* `EdgeCollapsedRips(arg; kwargs...)`: Use `arg` and `kwargs` to construct a `Rips`
  filtration and collapse it

* `EdgeCollapsedRips{I}(arg; kwargs...)`: Change the index type used to represent
  simplices. May be necessary for large inputs and high-dimensional computation.

# Examples

```jldoctest
julia> using Random; Random.seed!(1337);

julia> data = [tuple(rand(6)...) for _ in 1:100];

julia> rips = Rips(data)
Rips{Int64, Float64}(nv=100, sparse=false)

julia> length(TDARipserer.edges(rips))
3934

julia> collapsed = EdgeCollapsedRips(data) # or EdgeCollapsedRips(rips)
EdgeCollapsedRips{Int64, Float64}(nv=100)

julia> length(TDARipserer.edges(collapsed))
1324

julia> ripserer(rips) == ripserer(collapsed)
true

julia> ripserer(collapsed; dim_max=4)
5-element Vector{TDAPersistenceDiagrams.PersistenceDiagram}:
 100-element 0-dimensional PersistenceDiagram
 58-element 1-dimensional PersistenceDiagram
 35-element 2-dimensional PersistenceDiagram
 10-element 3-dimensional PersistenceDiagram
 4-element 4-dimensional PersistenceDiagram

```

# Reference

Boissonnat, J. D., & Pritam, S. (2019). [Edge Collapse and Persistence of Flag
Complexes](https://hal.inria.fr/hal-02395227/).
"""
struct EdgeCollapsedRips{I,T} <: AbstractRipsFiltration{I,T}
    adj::SparseMatrixCSC{T,Int}
    threshold::T
end

function EdgeCollapsedRips(rips::Rips{I,T}; verbose=false, threshold=nothing) where {I,T}
    adj = _core_graph(rips; verbose=verbose)
    isnothing(threshold) && (threshold = rips.threshold)
    return EdgeCollapsedRips{I,T}(adj, T(threshold))
end
function EdgeCollapsedRips(rips::EdgeCollapsedRips{I}; verbose=false, threshold=rips.threshold) where {I}
    return EdgeCollapsedRips{I}(rips; threshold=threshold)
end
function EdgeCollapsedRips{I}(dist_or_points; verbose=false, kwargs...) where {I}
    return EdgeCollapsedRips(Rips{I}(dist_or_points; kwargs...); verbose=verbose)
end
function EdgeCollapsedRips{I}(
    rips::EdgeCollapsedRips{<:Any,T}; threshold=rips.threshold
) where {I,T}
    return EdgeCollapsedRips{I,T}(rips.adj, T(threshold))
end
function EdgeCollapsedRips(arg; kwargs...)
    return EdgeCollapsedRips{Int}(arg; kwargs...)
end

threshold(rips::EdgeCollapsedRips) = rips.threshold
adjacency_matrix(rips::EdgeCollapsedRips) = rips.adj
