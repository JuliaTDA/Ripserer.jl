# Pointwise weights keep the diagram's filtration usable without an n-by-n matrix.
struct _EuclideanPointDistances{D} <: AbstractMatrix{Float64}
    points::Vector{SVector{D,Float64}}
end

Base.size(dist::_EuclideanPointDistances) = (length(dist.points), length(dist.points))

@inline function _point_distance(p, q)
    delta = p - q
    squared = sum(abs2, delta)
    return squared >= floatmin(Float64) ? sqrt(squared) : hypot(delta...)
end

@inline Base.getindex(dist::_EuclideanPointDistances, i::Int, j::Int) =
    _point_distance(dist.points[i], dist.points[j])

function _reduced_points(points, ::Val{D}) where {D}
    bound = sqrt(floatmax(Float64) / (8D))
    result = SVector{D,Float64}[]
    sizehint!(result, length(points))
    for point in points
        length(point) == D || throw(ArgumentError("points must have the same length"))
        all(x -> x isa Real && abs(x) <= bound && Float64(x) == x, point) ||
            throw(ArgumentError("reduced Rips requires finite, Float64-representable coordinates in the KD-tree distance range"))
        # Treat signed zeros as the same location before deduplicating.
        p = SVector{D,Float64}(Tuple(point))
        push!(result, map(x -> iszero(x) ? 0.0 : x, p))
    end
    if !allunique(result)
        @warn "points not unique"
        unique!(result)
    end
    return result
end

# Negated vertex IDs give the same equal-weight edge order as Simplex.
const _ReducedEdge = Tuple{Float64,Int,Int}

@inline _reduced_edge(points, u, v) =
    (_point_distance(points[u], points[v]), -max(u, v), -min(u, v))
@inline _edge_endpoints(edge::_ReducedEdge) = (-edge[2], -edge[3])

function _edge_simplex(::Type{I}, edge::_ReducedEdge) where {I}
    return Simplex{1}(I.(_edge_endpoints(edge)), edge[1])
end

function _find_lune!(lune, tree, points, edge::_ReducedEdge)
    u, v = _edge_endpoints(edge)
    empty!(lune)
    radius = _radius_search_bound(edge[1], length(points[u]))
    inrange!(lune, tree, points[u], radius)
    filter!(w -> _lune_witness(points, u, v, w, edge), lune)
    return lune
end

@inline function _lune_witness(points, u, v, w, edge::_ReducedEdge)
    return w != u && w != v &&
        _reduced_edge(points, u, w) < edge && _reduced_edge(points, v, w) < edge
end

function _rng_edges(tree, points, limit)
    rng_edges = _ReducedEdge[]
    # With no explicit threshold, compute the usual Rips radius during the scan.
    eccentricities = isnothing(limit) ? zeros(length(points)) : Float64[]
    neighbors, lune = Int[], Int[]
    for u in 2:length(points)
        nearby, _ = knn(tree, points[u], min(12, length(points)))
        last_witness = u
        candidates = if isnothing(limit)
            1:(u - 1)
        else
            empty!(neighbors)
            inrange!(neighbors, tree, points[u], _radius_search_bound(limit, length(points[u])))
        end
        for v in candidates
            v < u || continue
            edge = _reduced_edge(points, u, v)
            r = edge[1]
            if isnothing(limit)
                eccentricities[u] = max(eccentricities[u], r)
                eccentricities[v] = max(eccentricities[v], r)
            else
                r <= limit || continue
            end
            _lune_witness(points, u, v, last_witness, edge) && continue
            witness = findfirst(w -> _lune_witness(points, u, v, w, edge), nearby)
            if !isnothing(witness)
                last_witness = nearby[witness]
            else
                _find_lune!(lune, tree, points, edge)
                if isempty(lune)
                    push!(rng_edges, edge)
                else
                    last_witness = first(lune)
                end
            end
        end
    end
    radius = isnothing(limit) ? minimum(eccentricities) : limit
    filter!(edge -> edge[1] <= radius, rng_edges)
    return sort!(rng_edges), radius
end

function _push_next_edge!(heap, row, buffer, frontiers, candidates, tree, points, limit)
    # Double each row's search radius only when its current batch is exhausted.
    while isempty(buffer) && frontiers[row] < limit
        lower = frontiers[row]
        if iszero(lower)
            _, nearest = nn(tree, points[row], v -> v >= row)
            upper = max(2nearest, sqrt(floatmin(Float64)))
        else
            upper = 2lower
        end
        upper = min(upper, limit)
        frontiers[row] = upper
        empty!(candidates)
        inrange!(candidates, tree, points[row], _radius_search_bound(upper, length(points[row])))
        for v in candidates
            v < row || continue
            edge = _reduced_edge(points, row, v)
            lower < edge[1] <= upper && push!(buffer, edge)
        end
        sort!(buffer; rev=true)
    end
    isempty(buffer) || push!(heap, pop!(buffer))
    return nothing
end

function _lens_witness(tree, points, edge)
    edge[1] >= sqrt(floatmin(Float64)) || return 0
    u, v = _edge_endpoints(edge)
    w, _ = nn(tree, (points[u] + points[v]) / 2)
    _lune_witness(points, u, v, w, edge) || return 0
    # The angle test from Lemma 3.5, with slack for floating-point rounding.
    a = _point_distance(points[u], points[w]) / edge[1]
    b = _point_distance(points[v], points[w]) / edge[1]
    return (a^2 + b^2 - 1) / 2 < -0.87 * a * b - 32eps(Float64) ? w : 0
end

function _lune_components!(reps, lune, queue, tree, points, edge)
    empty!(reps)
    w = _lens_witness(tree, points, edge)
    if w != 0
        push!(reps, w)
        return reps
    end
    _find_lune!(lune, tree, points, edge)
    while !isempty(lune)
        push!(reps, pop!(lune))
        empty!(queue)
        push!(queue, last(reps))
        q = 1
        while q <= length(queue)
            w = queue[q]
            j = 1
            while j <= length(lune)
                z = lune[j]
                if _reduced_edge(points, w, z) < edge
                    push!(queue, z)
                    lune[j] = last(lune)
                    pop!(lune)
                else
                    j += 1
                end
            end
            q += 1
        end
    end
    return reps
end

function _xor_edges!(result, column, other)
    empty!(result)
    i = j = 1
    while i <= length(column) && j <= length(other)
        a, b = column[i], other[j]
        if a < b
            push!(result, a)
            i += 1
        elseif b < a
            push!(result, b)
            j += 1
        else
            i += 1
            j += 1
        end
    end
    for k in i:length(column)
        push!(result, column[k])
    end
    for k in j:length(other)
        push!(result, other[k])
    end
    return result
end

function _apparent_boundary(points, edge, w)
    u, v = _edge_endpoints(edge)
    a, b = minmax(_reduced_edge(points, u, w), _reduced_edge(points, v, w))
    return (a, b, edge)
end

function _reduced_h1(filtration::Rips{I}, tree, points, births, limit, cutoff) where {I}
    intervals = PersistenceInterval[]
    isempty(births) && return intervals
    buffers = [_ReducedEdge[] for _ in points]
    frontiers = zeros(length(points))
    heap = BinaryMinHeap{_ReducedEdge}()
    candidates, lune, reps, queue = Int[], Int[], Int[], Int[]
    # Non-tree RNG edges are the only births that can survive their own edge.
    rng_births = Set(births)
    paired = Dict{_ReducedEdge,Vector{_ReducedEdge}}()
    witnesses = Dict{_ReducedEdge,Int}()
    column, scratch = _ReducedEdge[], _ReducedEdge[]
    for row in 2:length(points)
        _push_next_edge!(heap, row, buffers[row], frontiers, candidates, tree, points, limit)
    end
    while !isempty(heap) && length(paired) < length(births)
        edge = pop!(heap)
        row, v = _edge_endpoints(edge)
        _lune_components!(reps, lune, queue, tree, points, edge)
        if length(reps) > 1
            apparent = _apparent_boundary(points, edge, first(reps))
            for i in 2:length(reps)
                w = reps[i]
                _xor_edges!(column, _apparent_boundary(points, edge, w), apparent)
                while !isempty(column)
                    pivot = last(column)
                    if haskey(paired, pivot)
                        _xor_edges!(scratch, column, paired[pivot])
                    elseif pivot in rng_births
                        paired[pivot] = copy(column)
                        if edge[1] - pivot[1] > cutoff
                            bs = _edge_simplex(I, pivot)
                            ds = Simplex{2}((I(row), I(v), I(w)), edge[1])
                            push!(intervals, PersistenceInterval(pivot[1], edge[1];
                                birth_simplex=bs, death_simplex=ds))
                        end
                        break
                    else
                        # A boundary's last edge cannot be a tree edge. Recreate
                        # the apparent column of a non-RNG edge only when needed.
                        x = get!(witnesses, pivot) do
                            lens = _lens_witness(tree, points, pivot)
                            lens != 0 && return lens
                            return last(_find_lune!(lune, tree, points, pivot))
                        end
                        _xor_edges!(scratch, column, _apparent_boundary(points, pivot, x))
                    end
                    column, scratch = scratch, column
                end
            end
        end
        length(paired) == length(births) && break
        _push_next_edge!(heap, row, buffers[row], frontiers, candidates, tree, points, limit)
    end
    for edge in births
        haskey(paired, edge) && continue
        Inf > cutoff || continue
        bs = _edge_simplex(I, edge)
        push!(intervals, PersistenceInterval(edge[1], Inf; birth_simplex=bs, death_simplex=nothing))
    end
    return intervals
end

_reduced_rips(::Type{Rips}, points::AbstractVector; kwargs...) =
    _reduced_rips(Rips{Int}, points; kwargs...)

function _reduced_rips(::Type{Rips{I}}, points::AbstractVector;
    alg, dim_max, field, reps, cutoff, verbose, implicit,
    metric=Euclidean(1e-12), threshold=nothing,
) where {I}
    field === Mod{2} || throw(ArgumentError("reduced Rips requires coefficients over F₂"))
    dim_max isa Integer && dim_max in (0, 1) ||
        throw(ArgumentError("reduced Rips supports dim_max=0 or 1"))
    any(d -> _reps(reps, d), 0:dim_max) &&
        throw(ArgumentError("reduced Rips supports barcodes without representatives"))
    cutoff >= 0 || throw(ArgumentError("reduced Rips requires a nonnegative cutoff"))
    metric isa Euclidean || throw(ArgumentError("reduced Rips requires Euclidean distances"))
    D = isempty(points) ? 2 : length(first(points))
    D in (2, 3) || throw(ArgumentError("reduced Rips requires 2D or 3D points"))
    if !isnothing(threshold)
        threshold >= 0 || throw(ArgumentError("reduced Rips requires a nonnegative threshold"))
        threshold = Float64(threshold)
    end
    return _reduced_rips(I, _reduced_points(points, Val(D)), threshold, dim_max, cutoff, verbose)
end

function _reduced_rips(::Type{I}, points::Vector{SVector{D,Float64}}, threshold,
    dim_max, cutoff, verbose,
) where {I,D}
    n = length(points)
    I <: Signed && n <= typemax(I) || throw(ArgumentError("invalid Rips simplex index type"))
    weights = _EuclideanPointDistances(points)
    # Clamp searches to a finite upper bound, retaining the requested diagram threshold.
    if n > 1
        lower = reduce((a, b) -> min.(a, b), points)
        upper = reduce((a, b) -> max.(a, b), points)
        diameter_bound = nextfloat(_point_distance(lower, upper))
        tree = KDTree(points)
        @prog_println verbose "Finding relative neighborhood graph..."
        rng_edges, limit = _rng_edges(tree, points, isnothing(threshold) ? nothing : min(threshold, diameter_bound))
    else
        rng_edges, limit = _ReducedEdge[], 0.0
    end
    thresh = isnothing(threshold) ? limit : threshold
    filtration = Rips{I,Float64,typeof(weights)}(weights, thresh)
    if n >= 2
        largest_dim = min(dim_max + 1, n - 1)
        index_overflow_check(simplex_type(filtration, largest_dim), Mod{2}, n,
            "reduced Rips simplex indices overflow. Try using a bigger index type.")
    end
    edges = [_edge_simplex(I, edge) for edge in rng_edges]
    h0, cycles, _ = zeroth_intervals(filtration, cutoff, verbose, Mod{2}, false, edges)
    dim_max == 0 && return [h0]
    intervals = if isempty(cycles)
        PersistenceInterval[]
    else
        @prog_println verbose "Computing reduced Rips H₁ intervals..."
        births = [(birth(edge), -Int(u), -Int(v)) for edge in cycles for (u, v) in (vertices(edge),)]
        _reduced_h1(filtration, tree, points, births, limit, cutoff)
    end
    h1 = PersistenceDiagram(sort!(intervals; by=persistence);
        threshold=thresh, dim=1, field=Mod{2}, filtration)
    return [h0, h1]
end

_reduced_rips(::Type, args...; kwargs...) =
    throw(ArgumentError("reduced Rips requires Euclidean point coordinates"))
