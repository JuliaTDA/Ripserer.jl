function _ripserer(
    ::Val{:redzed}, filtration::Rips, cutoff, verbose, field, dim_max, reps, implicit
)
    field === Mod{2} || throw(ArgumentError("RedZeD requires coefficients over F₂"))
    dim_max isa Integer && dim_max >= 0 || throw(ArgumentError("dim_max must be nonnegative"))
    cutoff >= 0 || throw(ArgumentError("RedZeD requires a nonnegative cutoff"))
    any(d -> _reps(reps, d), 0:dim_max) &&
        throw(ArgumentError("RedZeD provides simplex pairs only. Use cohomology for cocycles."))

    distances = adjacency_matrix(filtration)
    eltype(distances) <: Union{Float32,Float64} && !issparse(distances) ||
        throw(ArgumentError("RedZeD requires dense Float32 or Float64 Rips distances"))
    all(iszero, births(filtration)) ||
        throw(ArgumentError("RedZeD requires all vertex births to be zero"))
    n = nv(filtration)
    for i in 1:n, j in (i + 1):n
        d = distances[j, i]
        isfinite(d) && d > 0 || throw(ArgumentError("RedZeD distances must be finite and positive"))
    end
    isnan(threshold(filtration)) && throw(ArgumentError("RedZeD threshold must be ordered"))
    work_dim = Int(min(dim_max, max(0, n - 2)))
    # The reference implementation stores combinatorial indices in Int64 and
    # vertex labels in Int32, independently of the filtration's index type.
    k = min(n ÷ 2, work_dim + 2)
    n <= typemax(Int32) && binomial(big(n), k) <= typemax(Int64) ||
        throw(ArgumentError("RedZeD simplex index range exceeded"))

    @prog_println verbose "Computing RedZeD persistence..."
    result = PersistenceDiagram[]
    if n > 0 && threshold(filtration) >= 0
        bars, pairs = RedZeD.redzed(distances, work_dim;
            threshold=threshold(filtration), simplex_pairs=true)
        for d in 0:work_dim
            intervals = PersistenceInterval[]
            for ((b, e), (young, old)) in zip(bars[d], pairs[d])
                e - b > cutoff || continue
                bs = simplex(filtration, Val(d), Int.(young) .+ 1)
                ds = isempty(old) ? nothing : simplex(filtration, Val(d + 1), Int.(old) .+ 1)
                push!(intervals, PersistenceInterval(b, e; birth_simplex=bs, death_simplex=ds))
            end
            push!(result, _backend_diagram(filtration, field, d, intervals))
        end
    end
    for d in length(result):dim_max
        push!(result, _backend_diagram(filtration, field, d, PersistenceInterval[]))
    end
    return result
end

function _backend_diagram(filtration, field, dimension, intervals)
    return PersistenceDiagram(sort!(intervals; by=persistence);
        threshold=Float64(threshold(filtration)), dim=dimension, field, filtration)
end
