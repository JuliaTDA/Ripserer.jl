"""
    PersistencePairing

Finite Vietoris–Rips persistence pairing, represented by the critical edges
for each birth/death. `(0,0)` denotes a constant zero birth in H0. Fields
`birth_edges`, `death_edges`, `dim`, and `n_points` permit inspection and reuse.
Pairings are locally constant combinatorial data; recompute after optimization
steps. Holding them fixed across a pairing change does not compute the current PH.
"""
struct PersistencePairing
    birth_edges::Vector{Tuple{Int,Int}}
    death_edges::Vector{Tuple{Int,Int}}
    dim::Int
    n_points::Int
end

function _critical_edge(cell,D)
    vs = sort!(collect(vertices(cell)))
    length(vs)==1 && return (0,0)
    value = -Inf
    edge = (0,0)
    for j in 2:length(vs),i in 1:j-1
        a,b = vs[i],vs[j]
        if D[a,b]>value
            value,edge = D[a,b],(a,b)
        end
    end
    return edge
end

# Entirely combinatorial helper; AD extensions treat this function as inactive.
function _persistence_pairing(X,input::Symbol,dimension::Integer,modulus::Integer)
    dimension in (0,1) || throw(ArgumentError("differentiable Rips supports H0 and H1"))
    modulus==2 || throw(ArgumentError("differentiable Rips currently supports coefficients F2"))
    eltype(X)<:AbstractFloat || throw(ArgumentError("use floating-point input"))
    all(isfinite,X) || throw(ArgumentError("input must be finite"))
    if input==:points
        size(X,1)>=1 && size(X,2)>=1 || throw(ArgumentError("points must be a nonempty d×n matrix"))
        n = size(X,2)
        D = [norm(X[:,i]-X[:,j]) for i in 1:n,j in 1:n]
    elseif input==:distances
        n = size(X,1)
        n>=1 && size(X,2)==n || throw(DimensionMismatch("distance input must be nonempty and square"))
        all(>=(0),X) && all(iszero,diag(X)) ||
            throw(ArgumentError("distances must be nonnegative with a zero diagonal"))
        isapprox(X,transpose(X);atol=1e-12,rtol=1e-12) ||
            throw(ArgumentError("distances must be symmetric"))
        D = Matrix((X+transpose(X))/2)
    else
        throw(ArgumentError("input must be :points or :distances"))
    end
    n<=dimension+1 && return PersistencePairing(Tuple{Int,Int}[],Tuple{Int,Int}[],Int(dimension),n)
    diagram = ripserer(Rips,D;dim_max=dimension,modulus=modulus)[dimension+1]
    finite = filter(isfinite,diagram)
    births = [_critical_edge(birth_simplex(interval),D) for interval in finite]
    deaths = [_critical_edge(death_simplex(interval),D) for interval in finite]
    return PersistencePairing(births,deaths,Int(dimension),n)
end

"""
    persistence_pairing(X; input=:points, dim=1, modulus=2)

Compute finite Rips critical-edge pairs. Point coordinates are columns of a d×n
floating matrix. `input=:distances` accepts a dense symmetric nonnegative n×n
matrix with zero diagonal. The initial supported subset is H0/H1 over F2 with
Euclidean distances and no threshold/cutoff/collapse. Essential intervals are
excluded. The lexicographically first maximizing edge is selected at ties;
gradients there are branch subgradients, not a unique classical derivative.
"""
persistence_pairing(X::AbstractMatrix;input::Symbol=:points,dim::Integer=1,modulus::Integer=2) =
    _persistence_pairing(X,input,dim,modulus)

function _endpoint_value(X,edge,input)
    i,j = edge
    i==0 && return zero(eltype(X))
    input==:distances && return (X[i,j]+X[j,i])/2
    squared = zero(eltype(X))
    for k in axes(X,1)
        squared += (X[k,i]-X[k,j])^2
    end
    # Coincident points use zero subgradient rather than sqrt'(0)=Inf.
    return iszero(squared) ? zero(squared) : sqrt(squared)
end

function _pair_coordinates(X,pairs,input)
    values = Matrix{eltype(X)}(undef,length(pairs.birth_edges),2)
    for i in axes(values,1)
        values[i,1] = _endpoint_value(X,pairs.birth_edges[i],input)
        values[i,2] = _endpoint_value(X,pairs.death_edges[i],input)
    end
    return values
end

function _differentiable_persistence(X,input,dim,modulus)
    pairs = _persistence_pairing(X,input,dim,modulus)
    return _pair_coordinates(X,pairs,input)
end

"""
    differentiable_persistence(X; input=:points, dim=1, modulus=2)
    differentiable_persistence(X, pairs::PersistencePairing; input=:points)

Return an n_bars×2 matrix of finite diagram **births and deaths**, differentiable
with respect to point coordinates or dense distance filtration values. Supports
the precisely defined Rips subset in [`persistence_pairing`](@ref). Pairing changes
and ties are nonsmooth; away from them the gradient is the critical-edge derivative.
Empty diagrams return a 0×2 matrix and have zero loss/gradient.

The ChainRulesCore extension supports Zygote; the Enzyme extension marks only the
pairing computation inactive, leaving numeric endpoint extraction differentiable.
Optional AD extensions require Julia >=1.9 and the respective AD package loaded.
Distance gradients are split equally between symmetric matrix entries. Diagonal
entries are fixed zero. The two-argument form holds the supplied pairing fixed;
use it only for local calculations and recompute pairings between optimizer steps.

Reference: Carrière et al. (2021), *Optimizing persistent homology based functions*,
https://proceedings.mlr.press/v139/carriere21a.html.
"""
function differentiable_persistence(X::AbstractMatrix;
        input::Symbol=:points,dim::Integer=1,modulus::Integer=2)
    return _differentiable_persistence(X,input,dim,modulus)
end
function differentiable_persistence(X::AbstractMatrix,pairs::PersistencePairing;input::Symbol=:points)
    n = input==:points ? size(X,2) : size(X,1)
    n==pairs.n_points || throw(DimensionMismatch("point count differs from pairing"))
    input in (:points,:distances) || throw(ArgumentError("input must be :points or :distances"))
    return _pair_coordinates(X,pairs,input)
end

function _total_persistence(values,power)
    total = zero(eltype(values))
    for i in axes(values,1)
        total += (values[i,2]-values[i,1])^power
    end
    return total
end

"""
    topological_loss(X; input=:points, dim=1, modulus=2, power=2,
        target=nothing, negate=false)

Loss helper built from differentiable diagram endpoints. Without `target`, returns
the sum of finite persistences to `power` (negate to maximize/preserve features).
With a nonnegative scalar `target`, returns the squared discrepancy between that
sum and the target. `power` must be finite and >=1. Essential intervals are not
optimized. Combine this term with application-specific data/regularization losses;
maximizing persistence alone may merely expand geometric scale.

```julia
using TDARipserer, Zygote
X = [0.0 1.0 1.2 0.1; 0.0 0.0 1.1 1.0]
gradient = Zygote.gradient(x -> topological_loss(x;target=0.5),X)[1]
```
"""
function topological_loss(X::AbstractMatrix;input::Symbol=:points,dim::Integer=1,
        modulus::Integer=2,power::Real=2,target=nothing,negate::Bool=false)
    isfinite(power) && power>=1 || throw(ArgumentError("power must be finite and >=1"))
    target===nothing || (target isa Real && isfinite(target) && target>=0) ||
        throw(ArgumentError("target must be finite and nonnegative"))
    values = differentiable_persistence(X;input=input,dim=dim,modulus=modulus)
    total = _total_persistence(values,power)
    loss = target===nothing ? total : (total-target)^2
    return negate ? -loss : loss
end
