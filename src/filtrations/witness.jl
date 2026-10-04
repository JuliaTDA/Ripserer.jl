function _witness_distances(points, landmarks, metric)
    points = collect(points)
    ids = collect(Int, landmarks)
    isempty(points) && throw(ArgumentError("the witness cloud must be nonempty"))
    isempty(ids) && throw(ArgumentError("at least one landmark is required"))
    allunique(ids) && all(i -> 1 <= i <= length(points), ids) ||
        throw(ArgumentError("landmark indices must be distinct and in bounds"))
    d = [Float64(evaluate(metric, x, points[i])) for x in points, i in ids]
    return _check_witness_distances(d), ids
end

function _check_witness_distances(distances)
    d = Matrix{Float64}(distances)
    size(d, 1) > 0 && size(d, 2) > 0 ||
        throw(ArgumentError("witness-to-landmark distances must be nonempty"))
    all(x -> isfinite(x) && x >= 0, d) ||
        throw(ArgumentError("distances must be finite and nonnegative"))
    return d
end

function _witness_threshold(threshold)
    !isnan(threshold) && threshold >= 0 ||
        throw(ArgumentError("threshold must be nonnegative"))
    return Float64(threshold)
end

"""
    Witness(points, landmarks; metric=Euclidean(), max_dimension=2, threshold=Inf)
    Witness(witness_to_landmark_distances; max_dimension=2, threshold=Inf)

Relaxed weak witness filtration on selected landmarks (de Silva and Carlsson,
2004). `landmarks` contains distinct indices into `points`. In the matrix
form rows are witnesses and columns are landmarks. Vertex labels in the
filtration are landmark positions `1:length(landmarks)`; `flt.landmarks` maps
them back to the original point indices.

For a nonvertex simplex σ its weak-witness value is
`max(0, min_w(max_{l∈σ} d(w,l) - min_{l∉σ} d(w,l)))`.
Its filtration value is the maximum of that value and all facet values, so
every face is weakly witnessed at the same relaxation. All vertices start at
zero. If σ contains every landmark, the outside minimum is infinity.
`max_dimension` is the largest **simplex** dimension stored; computing H_d
requires dimension d+1 to include potential killing simplices.

Unlike a Rips/flag approximation, higher simplices have their own witnesses.
Construction is combinatorial in the landmark count and stores only simplices
below `threshold`; use [`LazyWitness`](@ref) for many landmarks. PH itself is
computed by the existing `ripserer` algorithms, including representatives.

Reference: https://doi.org/10.2312/SPBG/SPBG04/157-166.
"""
struct Witness{I,T} <: AbstractCustomFiltration{I,T}
    filtration::Custom{I,T}
    landmarks::Vector{Int}
    distances::Matrix{Float64}
end

function _witness(d, ids; max_dimension=2, threshold=Inf, verbose=false)
    max_dimension isa Integer && max_dimension >= 0 ||
        throw(ArgumentError("max_dimension must be a nonnegative integer"))
    threshold = _witness_threshold(threshold)
    n = size(d, 2)
    maxdim = min(max_dimension, n - 1)
    index_overflow_check(Tuple(reverse(collect((n - maxdim):n))))
    # AbstractCustomFiltration uses an edge adjacency even for a vertex-only complex.
    dicts = [Dict{Int,Float64}() for _ in 0:max(1, maxdim)]
    for i in 1:n
        dicts[1][i] = 0.0
    end
    for k in 2:(maxdim + 1)
        for subset in IterTools.subsets(1:n, k)
            vs = Tuple(reverse(subset))
            facet_values = [get(dicts[k - 1], index(Tuple(face)), Inf)
                            for face in IterTools.subsets(vs, k - 1)]
            lower = maximum(facet_values)
            lower <= threshold || continue
            inside = Set(vs)
            raw = Inf
            for w in axes(d, 1)
                outer = Inf
                inner = 0.0
                for l in 1:n
                    if l in inside
                        inner = max(inner, d[w, l])
                    else
                        outer = min(outer, d[w, l])
                    end
                end
                raw = min(raw, max(0.0, inner - outer))
            end
            value = max(lower, raw)
            value <= threshold && (dicts[k][index(vs)] = value)
        end
    end
    filtration = Custom{Int,Float64}(_adjacency_matrix(dicts), dicts, threshold)
    return Witness{Int,Float64}(filtration, ids, d)
end

function Witness(points::AbstractVector, landmarks; metric=Euclidean(), kwargs...)
    d, ids = _witness_distances(points, landmarks, metric)
    return _witness(d, ids; kwargs...)
end
function Witness(distances::AbstractMatrix; kwargs...)
    d = _check_witness_distances(distances)
    return _witness(d, collect(1:size(d, 2)); kwargs...)
end
simplex_dicts(w::Witness) = simplex_dicts(w.filtration)
adjacency_matrix(w::Witness) = adjacency_matrix(w.filtration)
threshold(w::Witness) = threshold(w.filtration)

"""
    LazyWitness(points, landmarks; nu=1, metric=Euclidean(), threshold=Inf)
    LazyWitness(witness_to_landmark_distances; nu=1, threshold=Inf)

Lazy witness flag filtration of de Silva and Carlsson (2004). For witness w,
let m_nu(w) be its `nu`-th nearest landmark distance, or zero when `nu=0`.
Edge (i,j) is born at `max(0, min_w(max(d(w,i),d(w,j))-m_nu(w)))`;
each higher simplex is born at its largest edge value. Vertices start at zero.
`nu` must lie in `0:number_of_landmarks`. Zero-valued edges are supported.

The matrix/landmark-label conventions match [`Witness`](@ref). Storage is
O(witnesses*landmarks + landmarks²), avoiding a full witness-cloud distance
matrix. Edge construction takes O(witnesses*landmarks²), and PH uses the
implicit flag-complex engine. This is an approximation to the full witness
complex, and its barcode need not equal a Rips barcode of the original cloud.

Reference: https://doi.org/10.2312/SPBG/SPBG04/157-166.
"""
struct LazyWitness{I,T} <: AbstractRipsFiltration{I,T}
    adj::Matrix{T}
    threshold::T
    landmarks::Vector{Int}
    nu::Int
end

function _lazy_witness(d, ids; nu=1, threshold=Inf, verbose=false)
    n = size(d, 2)
    nu isa Integer && 0 <= nu <= n ||
        throw(ArgumentError("nu must be an integer between zero and the landmark count"))
    threshold = _witness_threshold(threshold)
    offsets = nu == 0 ? zeros(size(d, 1)) :
        [partialsort(collect(view(d, w, :)), nu) for w in axes(d, 1)]
    adj = zeros(n, n)
    for j in 2:n, i in 1:(j - 1)
        value = Inf
        for w in axes(d, 1)
            value = min(value, max(0.0, max(d[w, i], d[w, j]) - offsets[w]))
        end
        adj[i, j] = adj[j, i] = value
    end
    return LazyWitness{Int,Float64}(adj, threshold, ids, nu)
end

function LazyWitness(points::AbstractVector, landmarks; metric=Euclidean(), kwargs...)
    d, ids = _witness_distances(points, landmarks, metric)
    return _lazy_witness(d, ids; kwargs...)
end
function LazyWitness(distances::AbstractMatrix; kwargs...)
    d = _check_witness_distances(distances)
    return _lazy_witness(d, collect(1:size(d, 2)); kwargs...)
end
adjacency_matrix(w::LazyWitness) = w.adj
threshold(w::LazyWitness) = w.threshold

# Rips treats a zero edge as an absent sparse edge. Lazy witness values can
# legitimately be zero, and its matrix is dense, so override those checks.
function unsafe_simplex(::Type{S}, w::LazyWitness, vs) where {S<:Simplex}
    dim(S) == 0 && return S(vs[1], 0.0)
    value = 0.0
    for i in eachindex(vs), j in (i + 1):length(vs)
        value = max(value, w.adj[vs[i], vs[j]])
    end
    return value <= threshold(w) ? S(index(vs), value) : nothing
end
function unsafe_cofacet(::Type{S}, w::LazyWitness, sx, vs, new_vertex) where {S<:Simplex}
    return unsafe_simplex(S, w, vs)
end
