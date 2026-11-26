# Run from the package root with julia --project=. --threads=1 benchmark/reduced_rips.jl
using TDARipserer, Random, Statistics, LinearAlgebra

BLAS.set_num_threads(1)

function endpoints(diagrams)
    return [sort!([(birth(i), death(i)) for i in d]) for d in diagrams]
end

function same_barcodes(a, b)
    return all(zip(a, b)) do (x, y)
        length(x) == length(y) || return false
        all(zip(x, y)) do ((b1, d1), (b2, d2))
            # The pointwise and matrix-based distance formulas round differently.
            isapprox(b1, b2; rtol=1e-8, atol=1e-12) &&
                (d1 == d2 || isapprox(d1, d2; rtol=1e-8, atol=1e-12))
        end
    end
end

function measure(work)
    result = endpoints(work())
    measurements = map(1:3) do _
        GC.gc()
        @timed work()
    end
    time = median(sample.time for sample in measurements)
    mib = median(sample.bytes for sample in measurements) / 2^20
    return result, round(time; sigdigits=4), round(mib; sigdigits=4)
end

function solid_torus(n)
    rng = MersenneTwister(691)
    points = NTuple{3,Float64}[]
    while length(points) < n
        x = 8rand(rng) - 4
        y = 8rand(rng) - 4
        z = 2rand(rng) - 1
        (hypot(x, y) - 3)^2 + z^2 <= 1 && push!(points, (x, y, z))
    end
    return points
end

function benchmark(label, points; threshold=nothing, collapse=false, graph=false)
    println(label, " (", length(points), " points)")
    baseline, seconds, mib = measure(() -> ripserer(points; threshold))
    println("  Rips                ", seconds, " s  ", mib, " MiB")
    methods = [
        "Rips sparse      " => () -> ripserer(points; threshold, sparse=true),
        "Reduced Rips     " => () -> ripserer(points; threshold, alg=:reduced_rips),
    ]
    if collapse
        push!(methods, "EdgeCollapsedRips" => () -> ripserer(EdgeCollapsedRips, points; threshold))
    end
    if graph
        push!(methods, "Rips radius graph" => () -> ripserer(points; threshold, radius_graph=true))
    end
    for (name, work) in methods
        result, seconds, mib = measure(work)
        @assert same_barcodes(baseline, result) name
        println("  ", name, "   ", seconds, " s  ", mib, " MiB")
        flush(stdout)
    end
end

println("Julia ", VERSION, ", ", Threads.nthreads(), " Julia thread(s); warm medians and cumulative allocations")
torus = solid_torus(5_000)
dragon = [Tuple(parse.(Float64, split(line)))::NTuple{3,Float64}
    for line in eachline(joinpath(@__DIR__, "data", "dragon2000.pts"))]
benchmark("Solid torus", torus[1:500]; collapse=true)
benchmark("Solid torus", torus[1:1_600])
benchmark("Dragon", dragon)
benchmark("Solid torus, threshold 0.4", torus; threshold=0.4, graph=true)
