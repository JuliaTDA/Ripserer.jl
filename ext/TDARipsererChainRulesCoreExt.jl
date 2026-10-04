module TDARipsererChainRulesCoreExt
using TDARipserer
using ChainRulesCore

function edge_gradient!(gradient,X,edge,weight,input)
    i,j = edge
    i==0 && return
    if input==:distances
        gradient[i,j] += weight/2
        gradient[j,i] += weight/2
    else
        distance = TDARipserer._endpoint_value(X,edge,input)
        iszero(distance) && return
        for k in axes(X,1)
            component = weight*(X[k,i]-X[k,j])/distance
            gradient[k,i] += component
            gradient[k,j] -= component
        end
    end
end

function coordinates_pullback(X,pairs,input,delta)
    delta = unthunk(delta)
    gradient = zeros(eltype(X),size(X))
    delta isa AbstractZero && return ProjectTo(X)(gradient)
    for i in eachindex(pairs.birth_edges)
        edge_gradient!(gradient,X,pairs.birth_edges[i],delta[i,1],input)
        edge_gradient!(gradient,X,pairs.death_edges[i],delta[i,2],input)
    end
    return ProjectTo(X)(gradient)
end

function ChainRulesCore.rrule(::typeof(TDARipserer._differentiable_persistence),X,input,dim,modulus)
    pairs = TDARipserer._persistence_pairing(X,input,dim,modulus)
    output = TDARipserer._pair_coordinates(X,pairs,input)
    pullback(delta) = (NoTangent(),coordinates_pullback(X,pairs,input,delta),
        NoTangent(),NoTangent(),NoTangent())
    return output,pullback
end

function ChainRulesCore.rrule(::typeof(TDARipserer._pair_coordinates),X,pairs,input)
    output = TDARipserer._pair_coordinates(X,pairs,input)
    pullback(delta) = (NoTangent(),coordinates_pullback(X,pairs,input,delta),
        NoTangent(),NoTangent())
    return output,pullback
end

function ChainRulesCore.rrule(::typeof(TDARipserer._total_persistence),values,power)
    output = TDARipserer._total_persistence(values,power)
    function pullback(delta)
        delta = unthunk(delta)
        gradient = zeros(eltype(values),size(values))
        if !(delta isa AbstractZero)
            for i in axes(values,1)
                weight = delta*power*(values[i,2]-values[i,1])^(power-1)
                gradient[i,1] = -weight
                gradient[i,2] = weight
            end
        end
        return NoTangent(),ProjectTo(values)(gradient),NoTangent()
    end
    return output,pullback
end
end
