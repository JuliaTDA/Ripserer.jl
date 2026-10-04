"""A certified optimal F₂ chain and a boundary certificate (`filling`) of homology."""
struct OptimalChainResult
    chain::Vector{Any}
    filling::Vector{Any}
    objective::Float64
    at::Float64
    dimension::Int
    tie_break::Symbol
end

_cellkey(s) = Tuple(sort!(collect(vertices(s))))

function _chain_cells(f::Union{Rips,Custom},d,at,max_cells)
    d>=0 && max_cells>0 || throw(ArgumentError("nonnegative dimension and positive max_cells required"))
    isfinite(at) || throw(ArgumentError("at must be finite"))
    at <= threshold(f) || throw(ArgumentError("at exceeds filtration threshold"))
    cells=[Any[] for _ in 0:d]
    n=nv(f); count=0
    function visit(v)
        sx=simplex(f,Val(length(v)-1),Tuple(reverse(v)))
        (isnothing(sx) || birth(sx)>at) && return
        count+=1
        count<=max_cells || throw(ArgumentError("max_cells exceeded"))
        push!(cells[length(v)],sx)
        length(v)==d+1 && return
        for j in last(v)+1:n
            visit(vcat(v,j))
        end
    end
    for i in 1:n
        visit([i])
    end
    return cells
end

function _mod2_boundary(rows,cols)
    lookup=Dict(_cellkey(s)=>i for (i,s) in enumerate(rows))
    I=Int[];J=Int[]
    for (j,s) in enumerate(cols)
        vs=_cellkey(s)
        for k in eachindex(vs)
            face=Tuple(vs[i] for i in eachindex(vs) if i!=k)
            haskey(lookup,face) || throw(ArgumentError("missing boundary face"))
            push!(I,lookup[face]);push!(J,j)
        end
    end
    return sparse(I,J,ones(Int,length(I)),length(rows),length(cols))
end

function _seed_vector(cells,chain)
    lookup=Dict(_cellkey(s)=>i for (i,s) in enumerate(cells))
    c=zeros(Int,length(cells))
    for s in chain
        key=_cellkey(s)
        haskey(lookup,key) || throw(ArgumentError("representative contains a cell absent at this threshold"))
        c[lookup[key]]=xor(c[lookup[key]],1)
    end
    return c
end

function _interval_chain(f,interval,at,kind)
    hasproperty(interval,:representative) || throw(ArgumentError("compute an interval with reps=true first"))
    reps=representative(interval)
    isempty(reps) && throw(ArgumentError("empty representative"))
    all(x -> coefficient(x) isa Mod{2},reps) || throw(ArgumentError("only F₂ representatives are supported"))
    if kind==:cocycle
        dim(simplex(first(reps)))==1 || throw(ArgumentError("cocycle conversion supports H₁ only; use homology reps with representative_kind=:cycle"))
        return Any[reconstruct_cycle(f,interval,at)...]
    elseif kind==:cycle
        return Any[simplex(x) for x in reps if !iszero(coefficient(x))]
    end
    throw(ArgumentError("representative_kind must be :cocycle or :cycle"))
end

function _weights(cells,weights)
    w = isnothing(weights) ? ones(length(cells)) : weights isa Function ? weights.(cells) : collect(weights)
    length(w)==length(cells) && all(x -> isfinite(x) && x>=0,w) || throw(ArgumentError("finite nonnegative weight per cell required"))
    return Float64.(w)
end

function _solve_mod2 end
_chain_optimizer_loaded() = isdefined(Base,:get_extension) &&
    Base.get_extension(@__MODULE__,:TDARipsererJuMPExt)!==nothing

"""
    optimal_cycle(f, interval; at=birth(interval), representative_kind=:cocycle,
                  weights=nothing, tie_break=:lexicographic, max_cells=10_000)

Minimum weighted support in the selected homology class over F₂, on Rips or
Custom simplicial filtrations. Requires `using JuMP, HiGHS`. The default converts
an H₁ cocycle using `reconstruct_cycle`; for homology representatives (including
higher dimensions), pass `representative_kind=:cycle`. `at` must be in the interval.
`weights` is a vector in enumerated cell order or a function of a simplex; default
unit weights minimize support size. The result's `filling` certifies z=c+∂y mod 2.
Lexicographic tie breaking chooses the smallest binary support vector among
solutions within `tolerance` of optimum, with one extra solve per cell. Use
`:solver` for any optimal support. Only a solver-certified optimum is returned.
"""
function optimal_cycle(f::Union{Rips,Custom},interval;at=birth(interval),
        representative_kind=:cocycle,weights=nothing,tie_break=:lexicographic,
        max_cells=10_000,tolerance=1e-8,kwargs...)
    birth(interval)<=at<death(interval) || throw(ArgumentError("at must lie in [birth,death)"))
    seed=_interval_chain(f,interval,at,representative_kind)
    d=dim(first(seed))
    cells=_chain_cells(f,d+1,at,max_cells)
    lower=d==0 ? Any[] : cells[d]
    target,higher=cells[d+1],cells[d+2]
    c=_seed_vector(target,seed)
    down=d==0 ? spzeros(Int,0,length(target)) : _mod2_boundary(lower,target)
    all(iseven,down*c) || throw(ArgumentError("representative is not a cycle"))
    B=_mod2_boundary(target,higher)
    _chain_optimizer_loaded() || throw(ArgumentError("chain optimization requires Julia 1.9+ with JuMP and HiGHS loaded"))
    z,y,cost=_solve_mod2(B,c,_weights(target,weights),:homologous;tie_break=tie_break,tolerance=tolerance,kwargs...)
    return OptimalChainResult(target[findall(z)],higher[findall(y)],cost,Float64(at),d,tie_break)
end

"""
    optimal_volume(f, cycle::AbstractVector; at, weights, tie_break, max_cells)
    optimal_volume(f, interval; at=death(interval), representative_kind=:cocycle, ...)

Minimum weighted F₂ filling whose boundary is the given cycle, via an integer
program. The interval form reconstructs its cycle at birth, then fills it at death.
An essential interval has no finite death and needs an explicit finite `at` with
a bounding chain. Infeasible fillings and non-optimal solver statuses are errors.
Supported filtrations, weights and tie policy are as for `optimal_cycle`.
"""
function optimal_volume(f::Union{Rips,Custom},cycle::AbstractVector;at,weights=nothing,
        tie_break=:lexicographic,max_cells=10_000,tolerance=1e-8,kwargs...)
    isempty(cycle) && throw(ArgumentError("cycle must be nonempty"))
    d=dim(first(cycle)); all(s -> dim(s)==d,cycle) || throw(ArgumentError("mixed-dimensional chain"))
    cells=_chain_cells(f,d+1,at,max_cells)
    target,higher=cells[d+1],cells[d+2]
    c=_seed_vector(target,cycle)
    down=d==0 ? spzeros(Int,0,length(target)) : _mod2_boundary(cells[d],target)
    all(iseven,down*c) || throw(ArgumentError("input is not a cycle"))
    _chain_optimizer_loaded() || throw(ArgumentError("chain optimization requires Julia 1.9+ with JuMP and HiGHS loaded"))
    z,_,cost=_solve_mod2(_mod2_boundary(target,higher),c,_weights(higher,weights),:filling;
        tie_break=tie_break,tolerance=tolerance,kwargs...)
    return OptimalChainResult(higher[findall(z)],Any[],cost,Float64(at),d+1,tie_break)
end
function optimal_volume(f::Union{Rips,Custom},interval;at=death(interval),representative_kind=:cocycle,kwargs...)
    return optimal_volume(f,_interval_chain(f,interval,birth(interval),representative_kind);at=at,kwargs...)
end
