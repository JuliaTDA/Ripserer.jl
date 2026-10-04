using TDARipserer, Flux, Zygote, Random, LinearAlgebra, Test

function flux_ph_example(;epochs=80,seed=20261003)
    rng = Xoshiro(seed)
    angles = collect(range(0,2π;length=13))[1:12]
    X = permutedims(hcat(cos.(angles),sin.(angles))) + 0.01randn(rng,2,12)
    model = Flux.Dense(0.4Matrix{Float64}(I,2,2),false)
    objective(m) = topological_loss(m(X);target=1.0,power=1)
    initial = objective(model)
    state = Flux.setup(Flux.Adam(0.02),model)
    for _ in 1:epochs
        value,grad = Flux.withgradient(objective,model)
        @assert isfinite(value)
        Flux.update!(state,model,grad[1])
    end
    final = objective(model)
    @test final < 0.1initial
    return (;initial,final,model)
end

if abspath(PROGRAM_FILE)==@__FILE__
    result = flux_ph_example()
    println("Flux topological loss: ",result.initial," → ",result.final)
end
