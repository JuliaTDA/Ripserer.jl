"""
    ZigzagEvent(operation, simplex)

A single `:add` or `:remove` of a simplex with distinct positive integer
vertex labels. Additions require all facets to be present; removals require
the simplex to be maximal. Vertex order is normalized. See
[`zigzag_persistence`](@ref) and [`ZigzagPersistence`](@ref).
"""
struct ZigzagEvent
    operation::Symbol
    simplex::_SimplexKey
    function ZigzagEvent(operation::Symbol, vertices)
        operation in (:add, :remove) || throw(ArgumentError("operation must be :add or :remove"))
        return new(operation, _simplex_key(vertices))
    end
end
ZigzagEvent(event::Pair) = ZigzagEvent(first(event), last(event))

const _ZigzagChanges = Vector{Pair{Int,Vector{_SimplexKey}}}

"""
    ZigzagRepresentative

An F₂ cycle along a zigzag bar, attached as `representative(interval)` when
`reps=true`. Use `cycle[state_index]` to retrieve the simplices with nonzero
coefficient at an observed state. `firstindex(cycle):lastindex(cycle)` is the
inclusive range of available state indices. Across every arrow of the bar,
the two cycles represent the same class after applying the inclusion.
Only changes in the chain are stored.
"""
struct ZigzagRepresentative
    first_state::Int
    last_state::Int
    changes::_ZigzagChanges
end

Base.firstindex(cycle::ZigzagRepresentative) = cycle.first_state
Base.lastindex(cycle::ZigzagRepresentative) = cycle.last_state

function _zz_chain_at(changes, t)
    i = searchsortedlast(changes, t => nothing; by=first)
    return iszero(i) ? _SimplexKey[] : last(changes[i])
end

function Base.getindex(cycle::ZigzagRepresentative, t::Integer)
    firstindex(cycle) <= t <= lastindex(cycle) || throw(BoundsError(cycle, t))
    return copy(_zz_chain_at(cycle.changes, t))
end

mutable struct _ZigzagSpace
    boundaries::Vector{BitSet}
    fillers::Vector{BitSet}
    cycles::Vector{BitSet}
    born::Vector{Int}
    basis::Union{Nothing,_BinaryBasis}
    paths::Union{Nothing,Vector{_ZigzagChanges}}
end
_ZigzagSpace(reps) = _ZigzagSpace(BitSet[], BitSet[], BitSet[], Int[], nothing,
    reps ? _ZigzagChanges[] : nothing)

"""
    ZigzagPersistence(; dim_max=1, initial=(), reps=false)

Online zigzag persistence over F2, with `push!(state, ZigzagEvent(...))`.
Sparse cycle, boundary and filling-chain bases are retained between updates;
no PH recomputation or all-pairs interval-rank calculation is performed.
The homology representatives are ordered by the right filtration of
Carlsson–de Silva: forward births append, backward births prepend, forward
deaths eliminate the last participating generator and backward deaths the
first. This works in arbitrary simplex/homology dimensions over F2.

`initial` specifies a complex by its maximal or all simplices; its face
closure is built at state index zero. Output includes dimensions `0:dim_max`.
The fields `steps`, `column_additions`, and `max_active_simplices` expose
reproducible work counters. Use [`zigzag_diagrams`](@ref) to inspect a stream
without closing it. Working storage depends on the active complex and its
cycle bases, not on the number of past boundary matrices.
Set `reps=true` to retain compatible cycle histories for the reported bars.
This mode retains chain changes from past events, including for completed bars.

Reference: https://arxiv.org/abs/0812.0197 and Carlsson, de Silva and Morozov,
*Zigzag persistent homology and real-valued functions* (SoCG 2009).
"""
mutable struct ZigzagPersistence
    active::Dict{_SimplexKey,Int}
    cofaces::Dict{Int,Set{Int}}
    spaces::Vector{_ZigzagSpace}
    intervals::Vector{Vector{PersistenceInterval}}
    dim_max::Int
    steps::Int
    next_id::Int
    free_ids::Vector{Int}
    column_additions::Int
    max_active_simplices::Int
    labels::Union{Nothing,Dict{Int,_SimplexKey}}
end

function _zz_space(state, d)
    while length(state.spaces) <= d
        push!(state.spaces, _ZigzagSpace(!isnothing(state.labels) && length(state.spaces) <= state.dim_max))
    end
    return state.spaces[d + 1]
end

function _zz_bar(space, i, death, last_state; kwargs...)
    b = space.born[i]
    isnothing(space.paths) && return PersistenceInterval(b, death; kwargs...)
    rep = ZigzagRepresentative(b, last_state, deepcopy(space.paths[i]))
    return PersistenceInterval(b, death; representative=rep, kwargs...)
end

function _zz_close!(state, d, space, i, t)
    if d <= state.dim_max && space.born[i] < t
        push!(state.intervals[d + 1], _zz_bar(space, i, t, t - 1))
    end
    !isnothing(space.paths) && deleteat!(space.paths, i)
    return nothing
end

function _zz_start_path!(state, space, t, chain; prepend=false)
    isnothing(space.paths) && return
    entry = _ZigzagChanges([t => sort!(_SimplexKey[state.labels[id] for id in chain])])
    if prepend
        pushfirst!(space.paths, entry)
    else
        push!(space.paths, entry)
    end
    return nothing
end

function _zz_add_path!(space, dest, src)
    isnothing(space.paths) && return
    a, b = space.paths[dest], space.paths[src]
    start = space.born[dest]
    times = sort!(union(first.(a), filter(>=(start), first.(b))))
    changes = _ZigzagChanges()
    for t in times
        chain = sort!(symdiff(_zz_chain_at(a, t), _zz_chain_at(b, t)))
        if isempty(changes) || chain != last(last(changes))
            push!(changes, t => chain)
        end
    end
    space.paths[dest] = changes
    return nothing
end

function _zz_add!(state, vs, t)
    haskey(state.active, vs) && throw(ArgumentError("simplex $vs is already present"))
    faces = _facets(vs)
    all(f -> haskey(state.active, f), faces) ||
        throw(ArgumentError("all facets must be present before adding $vs"))
    d = length(vs) - 1
    id = if isempty(state.free_ids)
        state.next_id += 1
        state.next_id - 1
    else
        pop!(state.free_ids)
    end
    state.active[vs] = id
    !isnothing(state.labels) && (state.labels[id] = vs)
    state.cofaces[id] = Set{Int}()
    for f in faces
        push!(state.cofaces[state.active[f]], id)
    end
    state.max_active_simplices = max(state.max_active_simplices, length(state.active))
    new_cycle = BitSet([id])
    if d > 0
        lower = _zz_space(state, d - 1)
        if isnothing(lower.basis)
            lower.basis = _binary_basis(vcat(lower.boundaries, lower.cycles))
            state.column_additions += lower.basis.additions
        end
        before = lower.basis.additions
        boundary = BitSet([state.active[f] for f in faces])
        remainder, coordinates = _binary_reduce(lower.basis, boundary)
        state.column_additions += lower.basis.additions - before
        isempty(remainder) || error("simplex boundary is outside the maintained cycle space")
        nb = length(lower.boundaries)
        homology = [i - nb for i in coordinates if i > nb]
        if !isempty(homology)
            # In a forward map, the last right-filtration summand participating
            # in the kernel dies; earlier summands retain their representatives.
            killed = maximum(homology)
            for i in homology
                i == killed || _zz_add_path!(lower, killed, i)
            end
            _zz_close!(state, d - 1, lower, killed, t)
            deleteat!(lower.cycles, killed)
            deleteat!(lower.born, killed)
            push!(lower.boundaries, boundary)
            push!(lower.fillers, BitSet([id]))
            lower.basis = nothing
            return
        end
        for i in coordinates
            symdiff!(new_cycle, lower.fillers[i])
            state.column_additions += 1
        end
    end
    space = _zz_space(state, d)
    push!(space.cycles, new_cycle)
    push!(space.born, t)
    _zz_start_path!(state, space, t, new_cycle)
    space.basis = nothing
end

function _zz_remove!(state, vs, t)
    haskey(state.active, vs) || throw(ArgumentError("simplex $vs is absent"))
    id = state.active[vs]
    isempty(state.cofaces[id]) || throw(ArgumentError("only maximal simplices can be removed"))
    d = length(vs) - 1
    space = _zz_space(state, d)
    killed = findfirst(c -> id in c, space.cycles)
    if !isnothing(killed)
        # The image of the backward inclusion consists of cycles whose
        # coefficient on the removed simplex is zero. Eliminate using the
        # earliest right-filtration representative that contains that simplex.
        pivot = space.cycles[killed]
        for i in eachindex(space.cycles)
            if i != killed && id in space.cycles[i]
                symdiff!(space.cycles[i], pivot)
                _zz_add_path!(space, i, killed)
                state.column_additions += 1
            end
        end
        # A filling chain may contain the removed simplex even when that
        # simplex also participates in homology. Add the pivot cycle to remove
        # its coefficient without changing the filling chain's boundary.
        if d > 0
            lower = _zz_space(state, d - 1)
            for filler in lower.fillers
                if id in filler
                    symdiff!(filler, pivot)
                    state.column_additions += 1
                end
            end
        end
        _zz_close!(state, d, space, killed, t)
        deleteat!(space.cycles, killed)
        deleteat!(space.born, killed)
        space.basis = nothing
    else
        d > 0 || error("a maximal vertex must represent a connected component")
        lower = _zz_space(state, d - 1)
        removed = findfirst(c -> id in c, lower.fillers)
        isnothing(removed) && error("missing filling chain for removed negative simplex")
        chain, cycle = lower.fillers[removed], lower.boundaries[removed]
        for i in eachindex(lower.fillers)
            if i != removed && id in lower.fillers[i]
                symdiff!(lower.fillers[i], chain)
                symdiff!(lower.boundaries[i], cycle)
                state.column_additions += 2
            end
        end
        deleteat!(lower.fillers, removed)
        deleteat!(lower.boundaries, removed)
        pushfirst!(lower.cycles, cycle)
        pushfirst!(lower.born, t)
        _zz_start_path!(state, lower, t, cycle; prepend=true)
        lower.basis = nothing
    end
    for f in _facets(vs)
        delete!(state.cofaces[state.active[f]], id)
    end
    delete!(state.cofaces, id)
    delete!(state.active, vs)
    !isnothing(state.labels) && delete!(state.labels, id)
    # All surviving representatives/fillers were cleared of this cell above.
    # Recycle chain coordinates so a long stream on a small complex does not
    # allocate ever-larger bitsets just because the event count grows.
    push!(state.free_ids, id)
end

function ZigzagPersistence(; dim_max=1, initial=(), reps::Bool=false)
    dim_max isa Integer && dim_max >= 0 || throw(ArgumentError("dim_max must be a nonnegative integer"))
    state = ZigzagPersistence(Dict{_SimplexKey,Int}(), Dict{Int,Set{Int}}(),
        _ZigzagSpace[], [PersistenceInterval[] for _ in 0:dim_max], dim_max, 0, 1, Int[], 0, 0,
        reps ? Dict{Int,_SimplexKey}() : nothing)
    for vs in _simplex_closure(initial)
        _zz_add!(state, vs, 0)
    end
    return state
end

function Base.push!(state::ZigzagPersistence, event::ZigzagEvent)
    t = state.steps + 1
    if event.operation == :add
        _zz_add!(state, event.simplex, t)
    else
        _zz_remove!(state, event.simplex, t)
    end
    state.steps = t
    return state
end
Base.push!(state::ZigzagPersistence, event::Pair) = push!(state, ZigzagEvent(event))

"""
    zigzag_diagrams(state::ZigzagPersistence)

Return a snapshot of the online barcode, as ordinary persistence diagrams.
State zero is the initial complex; state j follows event j. A finite interval
`[b,d)` is present at state indices b through d-1. Open intervals have death
`Inf` and are right-censored by the stream end, rather than asserting eternal
topological survival. Diagram metadata includes `kind=:zigzag` and `steps`.
When tracking representatives, each returned interval includes a snapshot
of its `ZigzagRepresentative`.
"""
function zigzag_diagrams(state::ZigzagPersistence)
    return map(0:state.dim_max) do d
        space = _zz_space(state, d)
        closed = state.intervals[d + 1]
        open = [_zz_bar(space, i, Inf, state.steps; right_censored=true) for i in eachindex(space.born)]
        return PersistenceDiagram(sort!(vcat(isnothing(space.paths) ? closed : deepcopy(closed), open));
            dim=d, kind=:zigzag, steps=state.steps)
    end
end

"""
    zigzag_persistence(events; dim_max=1, initial=(), reps=false)

Compute zigzag PH over F2 from simplex insertions/deletions. Each element is
a [`ZigzagEvent`](@ref) or `:add => vertices` / `:remove => vertices` pair.
Uses the same online sparse-basis engine as [`ZigzagPersistence`](@ref).

```julia
events = [:add => (1,), :add => (2,), :add => (1, 2), :remove => (1, 2)]
diagrams = zigzag_persistence(events; dim_max=0)
```
"""
function zigzag_persistence(events; dim_max=1, initial=(), reps::Bool=false)
    state = ZigzagPersistence(; dim_max=dim_max, initial=initial, reps=reps)
    for event in events
        push!(state, event)
    end
    return zigzag_diagrams(state)
end
