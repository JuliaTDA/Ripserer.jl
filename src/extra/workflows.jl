"""
    delay_embedding(series; dimension=3, delay=1)

DelayEmbeddings.embed adapter returning its native StateSpaceSet. Load
DelayEmbeddings first. Delays are in sample indices, not physical time; resample
irregular observations before use. See `suggest_delay` and the time-series example.
"""
function delay_embedding end

"""Suggest a delay using DelayEmbeddings.estimate_delay; this is guidance, not an optimum."""
function suggest_delay end

"""
    sliding_window_ph(points; window, step=1, dim_max=1, kwargs...)

Compute Rips PH of consecutive windows of *already embedded* states. Accepts a
vector of coordinate vectors or a matrix with states in columns. Returns a vector
of `(start, stop, diagrams)` records, where start/stop index embedded states.
Choose delay and dimension on training data; overlapping windows are dependent.
Remaining keywords are forwarded to ripserer. Window size must fit the embedding.
"""
function sliding_window_ph(points;window::Integer,step::Integer=1,dim_max::Integer=1,kwargs...)
    X=points isa AbstractMatrix ? collect(eachcol(points)) : collect(points)
    2<=window<=length(X) || throw(ArgumentError("window must be between 2 and the number of embedded states"))
    step>0 && dim_max>=0 || throw(ArgumentError("positive step and nonnegative dimension required"))
    d=length(first(X))
    all(x -> length(x)==d && all(isfinite,x),X) || throw(ArgumentError("finite, equally sized states required"))
    return [(start=i,stop=i+window-1,diagrams=ripserer(X[i:i+window-1];dim_max=dim_max,kwargs...))
        for i in 1:step:(length(X)-window+1)]
end

"""
    image_persistence(image; transform=identity, threshold=nothing, dim_max=1, kwargs...)

Apply a callable preprocessing transform (e.g. MetricSpaces.ImageFiltration), then
run Cubical PH. By default the threshold is the maximum *finite* transformed value,
excluding Inf background. NaN/-Inf and an empty finite foreground are rejected.
The image uses TDARipserer's vertex lower-star convention.
"""
function image_persistence(image;transform=identity,threshold=nothing,dim_max::Integer=1,kwargs...)
    A=transform(image)
    A isa AbstractArray{<:Real} && !isempty(A) || throw(ArgumentError("transform must return a nonempty real array"))
    all(x -> isfinite(x) || x==Inf,A) || throw(ArgumentError("NaN and -Inf filtration values are unsupported"))
    finite=filter(isfinite,vec(A))
    !isempty(finite) || throw(ArgumentError("transformed image has no finite vertices"))
    t=isnothing(threshold) ? maximum(finite) : threshold
    isfinite(t) || throw(ArgumentError("threshold must be finite"))
    return ripserer(Cubical(Array{Float64}(A);threshold=t);dim_max=dim_max,kwargs...)
end
