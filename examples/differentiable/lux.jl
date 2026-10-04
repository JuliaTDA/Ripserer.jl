using TDARipserer, Lux, Optimisers, Zygote, Random, LinearAlgebra, Test

function lux_ph_example(;epochs=80,seed=20261003)
    rng = Xoshiro(seed)
    angles = collect(range(0,2π;length=13))[1:12]
    X = permutedims(hcat(cos.(angles),sin.(angles))) + 0.01randn(rng,2,12)
    model = Lux.Dense(2=>2;use_bias=false)
    parameters,layer_state = Lux.setup(rng,model)
    parameters = merge(parameters,(weight=0.4Matrix{Float64}(I,2,2),))
    objective(p) = topological_loss(first(Lux.apply(model,X,p,layer_state));target=1.0,power=1)
    initial = objective(parameters)
    optimizer_state = Optimisers.setup(Optimisers.Adam(0.02),parameters)
    for _ in 1:epochs
        value,grad = Zygote.withgradient(objective,parameters)
        @assert isfinite(value)
        optimizer_state,parameters = Optimisers.update(optimizer_state,parameters,grad[1])
    end
    final = objective(parameters)
    @test final < 0.1initial
    return (;initial,final,model,parameters,layer_state)
end

if abspath(PROGRAM_FILE)==@__FILE__
    result = lux_ph_example()
    println("Lux topological loss: ",result.initial," → ",result.final)
end
