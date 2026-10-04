# Sparse exact linear algebra over F2. A chain is the set of indices with
# nonzero coefficient, so addition is symmetric difference.
const _SimplexKey = Tuple{Vararg{Int}}

function _simplex_key(vertices)
    vs = collect(Int, vertices)
    !isempty(vs) && all(>(0), vs) && allunique(vs) ||
        throw(ArgumentError("a simplex needs distinct positive integer vertices"))
    return Tuple(sort(vs))
end

_facets(vs::_SimplexKey) = length(vs) == 1 ? _SimplexKey[] :
    [_SimplexKey(Tuple(vs[j] for j in eachindex(vs) if j != i)) for i in eachindex(vs)]

mutable struct _BinaryBasis
    columns::Dict{Int,Tuple{BitSet,BitSet}}
    additions::Int
end
_BinaryBasis() = _BinaryBasis(Dict{Int,Tuple{BitSet,BitSet}}(), 0)

function _binary_reduce(basis::_BinaryBasis, chain::BitSet)
    remainder, coefficients = copy(chain), BitSet()
    while !isempty(remainder)
        pivot = maximum(remainder)
        haskey(basis.columns, pivot) || break
        column, tags = basis.columns[pivot]
        symdiff!(remainder, column)
        symdiff!(coefficients, tags)
        basis.additions += 1
    end
    return remainder, coefficients
end

function _binary_insert!(basis::_BinaryBasis, chain::BitSet, tag::Int)
    remainder, coefficients = _binary_reduce(basis, chain)
    isempty(remainder) && return false
    symdiff!(coefficients, BitSet([tag]))
    basis.columns[maximum(remainder)] = remainder, coefficients
    return true
end

function _binary_basis(chains)
    basis = _BinaryBasis()
    for (i, chain) in enumerate(chains)
        _binary_insert!(basis, chain, i) || error("dependent cycle basis")
    end
    return basis
end

function _simplex_closure(simplices)
    closed = Set{_SimplexKey}()
    function add(vs)
        vs in closed && return
        push!(closed, vs)
        foreach(add, _facets(vs))
    end
    for simplex in simplices
        add(_simplex_key(simplex))
    end
    return sort!(collect(closed); by=vs -> (length(vs), vs))
end
