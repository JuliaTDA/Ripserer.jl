# Independent dense F2 rank-invariant oracle. This intentionally does not use
# the production sparse cycle bases or interval decomposition algorithm.
function oracle_rref(matrix)
    a = Matrix{Bool}(matrix)
    pivots = Int[]
    row = 1
    for col in axes(a, 2)
        row > size(a, 1) && break
        pivot = findfirst(i -> a[i, col], row:size(a, 1))
        isnothing(pivot) && continue
        pivot += row - 1
        a[row, :], a[pivot, :] = copy(a[pivot, :]), copy(a[row, :])
        for i in axes(a, 1)
            i != row && a[i, col] && (a[i, :] .= xor.(a[i, :], a[row, :]))
        end
        push!(pivots, col)
        row += 1
    end
    return a, pivots
end
oracle_rank(a) = length(last(oracle_rref(a)))

function oracle_null(a)
    r, pivots = oracle_rref(a)
    free = setdiff(collect(axes(a, 2)), pivots)
    basis = zeros(Bool, size(a, 2), length(free))
    for (j, f) in enumerate(free)
        basis[f, j] = true
        for (row, p) in enumerate(pivots)
            basis[p, j] = r[row, f]
        end
    end
    return basis
end

function oracle_coordinates(basis, vector)
    r, pivots = oracle_rref(hcat(basis, vector))
    size(basis, 2) + 1 in pivots && error("oracle vector is not in the cycle span")
    answer = zeros(Bool, size(basis, 2))
    for (row, p) in enumerate(pivots)
        answer[p] = r[row, end]
    end
    return answer
end

function oracle_boundary(lower, upper)
    a = zeros(Bool, length(lower), length(upper))
    for (j, simplex) in enumerate(upper), (i, face) in enumerate(lower)
        length(face) + 1 == length(simplex) && all(v -> v in simplex, face) && (a[i, j] = true)
    end
    return a
end

function oracle_homology(complex, dimension)
    cells = sort([s for s in complex if length(s) == dimension + 1])
    lower = sort([s for s in complex if length(s) == dimension])
    upper = sort([s for s in complex if length(s) == dimension + 2])
    cycles = oracle_null(oracle_boundary(lower, cells))
    bd = oracle_boundary(cells, upper)
    boundaries = bd[:, last(oracle_rref(bd))]
    span, reps = boundaries, zeros(Bool, length(cells), 0)
    for j in axes(cycles, 2)
        extended = hcat(span, cycles[:, j])
        if oracle_rank(extended) > size(span, 2)
            span = extended
            reps = hcat(reps, cycles[:, j])
        end
    end
    return (; cells, span, reps, nb=size(boundaries, 2))
end

function oracle_inclusion(source, target)
    result = zeros(Bool, size(target.reps, 2), size(source.reps, 2))
    for j in axes(source.reps, 2)
        lifted = zeros(Bool, length(target.cells))
        for (i, simplex) in enumerate(source.cells)
            lifted[findfirst(==(simplex), target.cells)] = source.reps[i, j]
        end
        result[:, j] = oracle_coordinates(target.span, lifted)[(target.nb + 1):end]
    end
    return result
end

function oracle_generalized_rank(spaces, maps, directions, first_state, last_state)
    indices = first_state:last_state
    sizes = [size(spaces[i].reps, 2) for i in indices]
    offsets = cumsum(vcat(0, sizes))
    total = last(offsets)
    equations = zeros(Bool, 0, total)
    relations = zeros(Bool, total, 0)
    for i in first_state:(last_state - 1)
        local_index = i - first_state + 1
        source, target = directions[i] ? (local_index, local_index + 1) : (local_index + 1, local_index)
        a = maps[i]
        eq = zeros(Bool, sizes[target], total)
        rel = zeros(Bool, total, sizes[source])
        for j in 1:sizes[source]
            rel[offsets[source] + j, j] = true
        end
        for j in 1:sizes[target]
            eq[j, offsets[target] + j] = true
        end
        eq[:, (offsets[source] + 1):offsets[source + 1]] .= a
        rel[(offsets[target] + 1):offsets[target + 1], :] .= a
        equations = vcat(equations, eq)
        relations = hcat(relations, rel)
    end
    compatible = oracle_null(equations)
    image = zeros(Bool, total, size(compatible, 2))
    image[1:first(sizes), :] .= compatible[1:first(sizes), :]
    return oracle_rank(hcat(relations, image)) - oracle_rank(relations)
end
