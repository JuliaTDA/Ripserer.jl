"""
    VinePoint

One knot of a persistence vine: `time`, `birth`, `death`, and `dimension`.
Connect consecutive knots linearly. Infinite deaths stay infinite. Vines are
tracked through filtration-order transpositions, including zero-length bars.
"""
struct VinePoint
    time::Float64
    birth::Float64
    death::Float64
    dimension::Int
end

"""
    Vineyard(simplices; dim_max=1, time=0.0)
    Vineyard(filtration::AbstractFiltration; dim_max=1, time=0.0)

Incremental F2 persistence on a fixed simplicial complex. `simplices` is a
collection of `vertices => filtration_value` pairs. Missing faces are supplied
at the earliest coface value; explicitly specified face values must be <=
their coface values. The filtration constructor enumerates cells through
dimension `dim_max+1` using the existing filtration interface.

Retains a reduced boundary matrix R and upper triangular basis V satisfying
R = D*V. [`update_vineyard!`](@ref) follows affine filtration changes, updating
this factorization at adjacent order transpositions. It never reruns PH on a
new boundary matrix. `vines[id]` contains the knots of a continuous interval
track, and `swaps` records crossing times and work counts. `column_additions`
counts actual sparse F2 column additions; `transpositions` counts swaps.

The complex must stay fixed. Use [`ZigzagPersistence`](@ref) for cell
insertions/removals. At simultaneous equal-value crossings the endpoint
continuation is deterministic but need not be unique. No global matching of
independently recomputed diagrams is used.

Reference: Cohen-Steiner, Edelsbrunner and Morozov (2006),
https://mrzv.org/publications/vineyards/.
"""
mutable struct Vineyard
    simplices::Vector{_SimplexKey}
    values::Dict{_SimplexKey,Float64}
    reduced::Vector{BitSet}
    basis::Vector{BitSet}
    vine_ids::Dict{Tuple,Int}
    vines::Dict{Int,Vector{VinePoint}}
    swaps::Vector{NamedTuple}
    time::Float64
    dim_max::Int
    transpositions::Int
    column_additions::Int
end

function _vine_values(simplices)
    values = Dict{_SimplexKey,Float64}()
    for (vs, value) in simplices
        key = _simplex_key(vs)
        isfinite(value) || throw(ArgumentError("filtration values must be finite"))
        haskey(values, key) && values[key] != value &&
            throw(ArgumentError("conflicting values for simplex $key"))
        values[key] = Float64(value)
    end
    explicit = copy(values)
    function faces!(vs, value)
        for face in _facets(vs)
            if haskey(explicit, face)
                explicit[face] <= value || throw(ArgumentError("a face is born after its coface"))
            else
                values[face] = min(value, get(values, face, Inf))
            end
            faces!(face, values[face])
        end
    end
    for (vs, value) in explicit
        faces!(vs, value)
    end
    return values
end

function _reduce_binary_columns!(reduced, basis)
    owner = Dict{Int,Int}()
    additions = 0
    for j in eachindex(reduced)
        while !isempty(reduced[j])
            low = maximum(reduced[j])
            if haskey(owner, low)
                i = owner[low]
                symdiff!(reduced[j], reduced[i])
                symdiff!(basis[j], basis[i])
                additions += 1
            else
                owner[low] = j
                break
            end
        end
    end
    return additions
end

function _vine_pairs(state)
    deaths = Dict(maximum(c) => j for (j, c) in enumerate(state.reduced) if !isempty(c))
    return [(birth=state.simplices[i],
             death=haskey(deaths, i) ? state.simplices[deaths[i]] : nothing,
             dimension=length(state.simplices[i]) - 1)
            for i in eachindex(state.reduced) if isempty(state.reduced[i]) &&
                length(state.simplices[i]) - 1 <= state.dim_max]
end
_vine_key(pair) = (pair.birth, pair.death)

function _record_vines!(state, values, time; pairs=_vine_pairs(state),
                        ids=state.vine_ids, selected=nothing)
    for pair in pairs
        id = ids[_vine_key(pair)]
        isnothing(selected) || id in selected || continue
        knot = VinePoint(time, values[pair.birth],
            isnothing(pair.death) ? Inf : values[pair.death], pair.dimension)
        knots = state.vines[id]
        if isempty(knots) || (last(knots).time, last(knots).birth, last(knots).death) !=
                            (knot.time, knot.birth, knot.death)
            push!(knots, knot)
        end
    end
end

function Vineyard(simplices; dim_max=1, time=0.0)
    dim_max isa Integer && dim_max >= 0 || throw(ArgumentError("dim_max must be nonnegative"))
    isfinite(time) || throw(ArgumentError("time must be finite"))
    values = _vine_values(simplices)
    order = sort!(collect(keys(values)); by=s -> (values[s], length(s), s))
    positions = Dict(s => i for (i, s) in enumerate(order))
    reduced = [BitSet([positions[f] for f in _facets(s)]) for s in order]
    basis = [BitSet([i]) for i in eachindex(order)]
    additions = _reduce_binary_columns!(reduced, basis)
    state = Vineyard(order, values, reduced, basis, Dict{Tuple,Int}(),
        Dict{Int,Vector{VinePoint}}(), NamedTuple[], Float64(time), dim_max, 0, additions)
    for (id, pair) in enumerate(_vine_pairs(state))
        state.vine_ids[_vine_key(pair)] = id
        state.vines[id] = VinePoint[]
    end
    _record_vines!(state, values, Float64(time))
    return state
end

function Vineyard(flt::AbstractFiltration; dim_max=1, time=0.0)
    simplex_type(flt, 0) <: AbstractSimplex ||
        throw(ArgumentError("Vineyard requires a simplicial filtration; cubical cells are not supported"))
    pairs = Pair[(i,) => b for (i, b) in enumerate(births(flt))]
    cells = collect(edges(flt))
    append!(pairs, [Tuple(vertices(s)) => birth(s) for s in cells])
    for _ in 2:(dim_max + 1)
        isempty(cells) && break
        cells = unique(collect(columns_to_reduce(flt, cells)))
        append!(pairs, [Tuple(vertices(s)) => birth(s) for s in cells])
    end
    return Vineyard(pairs; dim_max=dim_max, time=time)
end

function _swap_rows!(column, i)
    a, b = i in column, i + 1 in column
    if a != b
        delete!(column, a ? i : i + 1)
        push!(column, a ? i + 1 : i)
    end
end

function _vine_swap!(state, i)
    left, right = state.simplices[i], state.simplices[i + 1]
    oldpairs, oldids = _vine_pairs(state), state.vine_ids
    additions = 0
    # Remove the only entry that would become subdiagonal after conjugation.
    if i in state.basis[i + 1]
        symdiff!(state.basis[i + 1], state.basis[i])
        symdiff!(state.reduced[i + 1], state.reduced[i])
        additions += 1
    end
    for column in state.reduced
        _swap_rows!(column, i)
    end
    for column in state.basis
        _swap_rows!(column, i)
    end
    state.reduced[i], state.reduced[i + 1] = state.reduced[i + 1], state.reduced[i]
    state.basis[i], state.basis[i + 1] = state.basis[i + 1], state.basis[i]
    state.simplices[i], state.simplices[i + 1] = right, left
    # Only reduced columns affected by the transposition need additions. The
    # owner scan detects those collisions, retaining every unaffected column.
    additions += _reduce_binary_columns!(state.reduced, state.basis)
    newids = Dict{Tuple,Int}()
    available = Set(values(oldids))
    newpairs = _vine_pairs(state)
    for pair in newpairs
        key = _vine_key(pair)
        if haskey(oldids, key)
            newids[key] = oldids[key]
            delete!(available, oldids[key])
        end
    end
    swapped_dimension = length(left) - 1
    for pair in newpairs
        key = _vine_key(pair)
        haskey(newids, key) && continue
        # At a same-dimensional cell swap, H_d birth endpoints cross and
        # H_(d-1) death endpoints cross. Continue along the other endpoint.
        candidates = [old for old in oldpairs if oldids[_vine_key(old)] in available &&
            old.dimension == pair.dimension &&
            (pair.dimension == swapped_dimension ? old.death == pair.death : old.birth == pair.birth)]
        isempty(candidates) && error("cannot continue vine through adjacent transposition")
        chosen = first(candidates)
        newids[key] = oldids[_vine_key(chosen)]
        delete!(available, newids[key])
    end
    isempty(available) || error("vine count changed on a fixed complex")
    state.vine_ids = newids
    state.transpositions += 1
    state.column_additions += additions
    changed = Set(oldids[key] for key in keys(oldids)
                  if !haskey(newids, key) || newids[key] != oldids[key])
    return (; additions, previous_pairs=oldpairs, previous_ids=oldids, changed)
end

"""
    update_vineyard!(state, values; time=state.time+1)

Update a [`Vineyard`](@ref) along the affine interpolation from its previous
filtration to `values`. Supply a dictionary/collection of simplex-value pairs
(partial updates are allowed), or a vector of values in the current
`state.simplices` order. The time must be finite and strictly increasing.
Validate all face inequalities before mutating the state.

Crossings are processed chronologically by adjacent transpositions. The
factorization is conjugated and repaired by sparse column additions; the
unaffected reduced columns are reused. Runtime depends on the number of
order crossings, so a large permutation can cost more than fresh PH. Updates
that preserve the order perform no reduction work. With coincident crossing
times, the leftmost available transposition is selected first.
"""
function update_vineyard!(state::Vineyard, updates; time=state.time + 1)
    isfinite(time) && time > state.time || throw(ArgumentError("time must be finite and increasing"))
    target = copy(state.values)
    entries = if updates isa AbstractVector{<:Real}
        length(updates) == length(state.simplices) || throw(DimensionMismatch("one value per simplex is required"))
        zip(state.simplices, updates)
    else
        updates
    end
    for (simplex, value) in entries
        key = _simplex_key(simplex)
        haskey(target, key) || throw(ArgumentError("updates cannot change the complex"))
        isfinite(value) || throw(ArgumentError("filtration values must be finite"))
        target[key] = Float64(value)
    end
    for (simplex, value) in target, face in _facets(simplex)
        target[face] <= value || throw(ArgumentError("a face is born after its coface"))
    end
    final_order = sort!(collect(keys(target)); by=s -> (target[s], length(s), s))
    rank = Dict(s => i for (i, s) in enumerate(final_order))
    start = copy(state.values)
    alpha = 0.0
    while state.simplices != final_order
        next_alpha, next_index = Inf, 0
        for i in 1:(length(state.simplices) - 1)
            a, b = state.simplices[i], state.simplices[i + 1]
            rank[a] > rank[b] || continue
            slope = (target[a] - start[a]) - (target[b] - start[b])
            crossing = iszero(slope) ? alpha : (start[b] - start[a]) / slope
            crossing = clamp(crossing, alpha, 1.0)
            if crossing < next_alpha
                next_alpha, next_index = crossing, i
            end
        end
        next_index > 0 || error("no adjacent inversion in a changed filtration order")
        alpha = next_alpha
        crossing_time = state.time + alpha * (time - state.time)
        crossing_values = Dict(s => start[s] + alpha * (target[s] - start[s]) for s in keys(start))
        a, b = state.simplices[next_index], state.simplices[next_index + 1]
        result = _vine_swap!(state, next_index)
        _record_vines!(state, crossing_values, crossing_time;
            pairs=result.previous_pairs, ids=result.previous_ids, selected=result.changed)
        _record_vines!(state, crossing_values, crossing_time; selected=result.changed)
        push!(state.swaps, (; time=crossing_time, left=a, right=b, column_additions=result.additions))
    end
    state.values = target
    state.time = Float64(time)
    _record_vines!(state, target, state.time)
    return state
end

"""
    vineyard_diagrams(state; include_zero=false)

Current persistence diagrams of a [`Vineyard`](@ref), with `vine_id`,
`birth_simplex` and `death_simplex` metadata on each interval. Zero-length
intervals are omitted by default, as in `ripserer`; they remain tracked in
`state.vines` so diagonal crossings do not erase an interval's identity.
"""
function vineyard_diagrams(state::Vineyard; include_zero=false)
    intervals = [PersistenceInterval[] for _ in 0:state.dim_max]
    for pair in _vine_pairs(state)
        b = state.values[pair.birth]
        d = isnothing(pair.death) ? Inf : state.values[pair.death]
        (include_zero || b < d) || continue
        push!(intervals[pair.dimension + 1], PersistenceInterval(b, d;
            vine_id=state.vine_ids[_vine_key(pair)], birth_simplex=pair.birth, death_simplex=pair.death))
    end
    return [PersistenceDiagram(sort(intervals[d + 1]); dim=d, kind=:vineyard, time=state.time)
            for d in 0:state.dim_max]
end
