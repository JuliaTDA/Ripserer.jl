function _cubical_h0(filtration::Cubical{N,T}, cutoff) where {N,T}
    data = filtration.data
    grid = CartesianIndices(data)
    linear = LinearIndices(data)
    direction = ntuple(i -> _one_hot(i, Val(N)), N)
    one_vertex = CartesianIndex(ntuple(_ -> 1, N))
    order = sortperm(vec(data))
    components = IntDisjointSets(length(data))
    elder = collect(eachindex(order))
    present = falses(length(data))
    intervals = PersistenceInterval[]
    for v in order
        level = data[v]
        level <= threshold(filtration) || break
        present[v] = true
        x = grid[v]
        for axis in direction, step in (-1, 1)
            y = x + step * axis
            checkbounds(Bool, data, y) || continue
            u = linear[y]
            present[u] || continue
            a, b = find_root!(components, v), find_root!(components, u)
            a == b && continue
            old, young = elder[a], elder[b]
            if (data[young], young) < (data[old], old)
                old, young = young, old
            end
            parent = root_union!(components, a, b)
            elder[parent] = old
            if Float64(level) - Float64(data[young]) > cutoff
                bs = Cube{0,T,N}(2 * grid[young] - one_vertex, data[young])
                ds = Cube{1,T,N}(x + y - one_vertex, level)
                push!(intervals, PersistenceInterval(birth(bs), level; birth_simplex=bs, death_simplex=ds))
            end
        end
    end
    for u in order
        present[u] || break
        find_root!(components, u) == u || continue
        v = elder[u]
        b = data[v]
        Inf - Float64(b) > cutoff || continue
        bs = Cube{0,T,N}(2 * grid[v] - one_vertex, b)
        push!(intervals, PersistenceInterval(b, Inf; birth_simplex=bs, death_simplex=nothing))
    end
    return intervals
end

struct _DualFace{N,T}
    first::Int
    second::Int
    root::CartesianIndex{N}
    grade::T
end

function _join_dual!(components, elder, grades, exterior, u, v)
    a, b = find_root!(components, u), find_root!(components, v)
    a == b && return 0
    old, young = elder[a], elder[b]
    if _older_dual(grades, exterior, young, old)
        old, young = young, old
    end
    elder[root_union!(components, a, b)] = old
    return young
end

function _older_dual(grades, exterior, u, v)
    return u == exterior || (v != exterior && (grades[u], u) > (grades[v], v))
end

function _cubical_dual(filtration::Cubical{N,T}, cutoff) where {N,T}
    intervals = PersistenceInterval[]
    sizes = size(filtration.data) .- 1
    any(<=(0), sizes) && return intervals
    top = CartesianIndices(sizes)
    linear = LinearIndices(top)
    direction = ntuple(i -> _one_hot(i, Val(N)), N)
    grades = [filtration.cubemap[2 * x] for x in top]
    exterior = length(top) + 1
    components = IntDisjointSets(exterior)
    elder = collect(1:exterior)

    # Pair a top cell with an incident equal-grade face leading to an older
    # cell or to the exterior. Every pair is a zero bar and the strict age
    # order makes the chosen pairs a forest. Collapse it before sorting faces.
    for x in top
        u = linear[x]
        for axis in direction, step in (-1, 1)
            y = x + step * axis
            v = checkbounds(Bool, top, y) ? linear[y] : exterior
            face = 2 * x + step * axis
            if filtration.cubemap[face] == grades[u] && _older_dual(grades, exterior, v, u)
                _join_dual!(components, elder, grades, exterior, u, v)
                break
            end
        end
    end

    faces = _DualFace{N,T}[]
    for x in top, axis in direction
        u = linear[x]
        y = x + axis
        v = checkbounds(Bool, top, y) ? linear[y] : exterior
        if find_root!(components, u) != find_root!(components, v)
            index = 2 * x + axis
            push!(faces, _DualFace(u, v, index, filtration.cubemap[index]))
        end
        if !checkbounds(Bool, top, x - axis) &&
            find_root!(components, u) != find_root!(components, exterior)
            index = 2 * x - axis
            push!(faces, _DualFace(u, exterior, index, filtration.cubemap[index]))
        end
    end
    sort!(faces; by=face -> face.grade, rev=true)
    for face in faces
        u = _join_dual!(components, elder, grades, exterior, face.first, face.second)
        b = face.grade
        (u == 0 || b > threshold(filtration)) && continue
        d = grades[u]
        censored = d > threshold(filtration)
        e = censored ? Inf : d
        Float64(e) - Float64(b) > cutoff || continue
        bs = Cube{N - 1,T,N}(face.root, b)
        ds = censored ? nothing : Cube{N,T,N}(2 * top[u], d)
        push!(intervals, PersistenceInterval(b, e; birth_simplex=bs, death_simplex=ds))
    end
    return intervals
end

function _ripserer(
    ::Val{:cubical_union_find}, filtration::Cubical{N}, cutoff, verbose, field, dim_max,
    reps, implicit
) where {N}
    if !(N in (2, 3)) || field !== Mod{2} || cutoff < 0 ||
        any(d -> _reps(reps, d), 0:dim_max) || !all(isfinite, filtration.data) ||
        isnan(threshold(filtration))
        return _ripserer(Val(:cohomology), filtration, cutoff, verbose, field, dim_max, reps, implicit)
    end
    result = if N == 3 && dim_max >= 1
        _ripserer(Val(:cohomology), filtration, cutoff, verbose, field, 1, reps, implicit)
    else
        [_backend_diagram(filtration, field, 0, _cubical_h0(filtration, cutoff))]
    end
    if dim_max >= N - 1
        @prog_println verbose "Computing dual cubical persistence..."
        push!(result, _backend_diagram(filtration, field, N - 1, _cubical_dual(filtration, cutoff)))
    end
    for d in length(result):dim_max
        push!(result, _backend_diagram(filtration, field, d, PersistenceInterval[]))
    end
    return result
end
