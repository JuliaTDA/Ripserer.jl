# Run: julia --project=TDARipserer.jl TDARipserer.jl/benchmark/p3.jl
# Native Julia timing/allocation tools keep this runnable without extra deps.
using TDARipserer
using TDAPersistenceDiagrams

function measure(label, work; repetitions=3)
    work() # compile and warm the same workload
    times, bytes = Float64[], Int[]
    for _ in 1:repetitions
        GC.gc()
        push!(bytes, @allocated work())
        push!(times, @elapsed work())
    end
    println(label, ": best_seconds=", minimum(times), " min_allocated_bytes=", minimum(bytes))
end

circle(n) = [[cos(t), sin(t)] for t in range(0, 2π; length=n + 1)[1:n]]

function witness_benchmark()
    points = circle(1200)
    landmarks = collect(1:25:1200)
    witness_work = () -> ripserer(LazyWitness(points, landmarks; threshold=0.4))
    rips_work = () -> ripserer(Rips(points; threshold=0.4))
    @assert !isempty(witness_work()[2]) && !isempty(rips_work()[2])
    measure("1200 points: lazy witness with 48 landmarks", witness_work)
    measure("1200 points: full Rips", rips_work)

    # Large workload whose dense full-cloud matrix alone would require 3.2 GB.
    large = circle(20_000)
    landmarks = collect(round.(Int, range(1, length(large); length=64)))
    measure("20000 witnesses / 64 landmarks", () -> begin
        result = ripserer(LazyWitness(large, landmarks; threshold=0.4))
        @assert !isempty(result[2])
    end; repetitions=2)
    println("large dense Rips distance matrix bytes avoided=", 8length(large)^2,
            " witness-to-landmark matrix bytes=", 8length(large) * length(landmarks))
end

function zigzag_benchmark()
    events = [ZigzagEvent(isodd(i) ? :add : :remove, (1, 2, 3)) for i in 1:10_000]
    function run()
        state = ZigzagPersistence(; initial=[(1, 2), (2, 3), (1, 3)])
        for event in events
            push!(state, event)
        end
        return state
    end
    measure("10000 zigzag triangle additions/deletions", run)
    state = run()
    println("zigzag column_additions=", state.column_additions,
            " max_active_simplices=", state.max_active_simplices)
end

function vineyard_benchmark()
    n = 60
    simplices = vcat([(i,) => 0.0 for i in 1:n],
        [(i, mod1(i + 1, n)) => (2 + sin(i)) for i in 1:n])
    frames = [Dict(Tuple(sort(collect(first(p)))) =>
        (length(first(p)) == 1 ? 0.0 : 2 + sin(first(first(p)) + frame / 20))
        for p in simplices) for frame in 1:30]
    function incremental()
        state = Vineyard(simplices)
        for (time, frame) in enumerate(frames)
            update_vineyard!(state, frame; time=time)
        end
        return state
    end
    measure("30 vineyard frames / 120 cells", incremental)
    measure("30 fresh reductions / 120 cells", () -> [Vineyard(collect(f)) for f in frames])
    state = incremental()
    println("vineyard transpositions=", state.transpositions,
            " column_additions=", state.column_additions)
end

println("Julia ", VERSION, "; deterministic circle and moving-cycle workloads")
witness_benchmark()
zigzag_benchmark()
vineyard_benchmark()
