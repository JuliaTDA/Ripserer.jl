module TDARipsererDelayEmbeddingsExt
using TDARipserer
using DelayEmbeddings

function TDARipserer.delay_embedding(series::AbstractVector{<:Real};dimension::Integer=3,delay::Integer=1)
    dimension>0 && delay>0 || throw(ArgumentError("positive embedding dimension and delay required"))
    all(isfinite,series) || throw(ArgumentError("series must be finite"))
    length(series)>(dimension-1)*delay || throw(ArgumentError("series too short for embedding"))
    return DelayEmbeddings.embed(series,dimension,delay)
end
function TDARipserer.suggest_delay(series::AbstractVector{<:Real};method="mi_min",lags=1:100,kwargs...)
    length(series)>2 && all(isfinite,series) || throw(ArgumentError("finite nontrivial series required"))
    maximum(series)>minimum(series) || throw(ArgumentError("constant series has no identifiable delay"))
    all(x -> x isa Integer && x>0,lags) && !isempty(lags) || throw(ArgumentError("positive lag indices required"))
    maximum(lags)<length(series) || throw(ArgumentError("lags exceed series length"))
    return DelayEmbeddings.estimate_delay(series,method,lags;kwargs...)
end
end
