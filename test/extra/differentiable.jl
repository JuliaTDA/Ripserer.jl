using TDARipserer, TDAPersistenceDiagrams, Test, LinearAlgebra, Random, Zygote

function central_gradient(f,X;epsilon=1e-6)
    gradient = similar(X)
    for i in eachindex(X)
        plus,minus = copy(X),copy(X)
        plus[i]+=epsilon
        minus[i]-=epsilon
        gradient[i]=(f(plus)-f(minus))/(2epsilon)
    end
    return gradient
end

@testset "finite Rips endpoint gradients" begin
    rng = Xoshiro(20261003)
    X = [0.0 1.0 1.23 0.11; 0.0 0.07 1.1 1.02]
    for dim in (0,1)
        pairs = persistence_pairing(X;dim=dim)
        coordinates = differentiable_persistence(X;dim=dim)
        diagram = filter(isfinite,ripserer([Tuple(X[:,i]) for i in axes(X,2)];dim_max=dim)[dim+1])
        @test coordinates[:,1] ≈ birth.(diagram)
        @test coordinates[:,2] ≈ death.(diagram)
        @test differentiable_persistence(X,pairs) ≈ coordinates
        weights = randn(rng,size(coordinates))
        f = x -> sum(weights .* differentiable_persistence(x;dim=dim))
        analytic = Zygote.gradient(f,X)[1]
        numerical = central_gradient(f,X)
        @test analytic ≈ numerical atol=1e-6 rtol=1e-5
        # Birth and death must propagate independently, not only their difference.
        for column in (1,2)
            endpoint = x -> sum(differentiable_persistence(x;dim=dim)[:,column])
            @test Zygote.gradient(endpoint,X)[1] ≈ central_gradient(endpoint,X) atol=1e-6
        end
        @test vec(sum(analytic;dims=2)) ≈ zeros(2) atol=1e-10 # translation invariance
    end
    loss = x -> topological_loss(x;dim=1,target=0.7,power=2)
    @test Zygote.gradient(loss,X)[1] ≈ central_gradient(loss,X) atol=1e-6
    @test topological_loss(X;negate=true)==-topological_loss(X)
    @test topological_loss(Float32.(X)) isa Float32

    D = [norm(X[:,i]-X[:,j]) for i in axes(X,2),j in axes(X,2)]
    f = d -> topological_loss(d;input=:distances,dim=1)
    gradient = Zygote.gradient(f,D)[1]
    @test gradient ≈ transpose(gradient)
    direction = randn(rng,size(D)); direction = (direction+transpose(direction))/2
    direction[diagind(direction)].=0
    @test sum(gradient .* direction) ≈ (f(D+1e-6direction)-f(D-1e-6direction))/2e-6 atol=1e-6
    empty = X[:,1:2]
    @test size(differentiable_persistence(empty))==(0,2)
    @test topological_loss(empty)==0
    @test all(iszero,Zygote.gradient(x -> topological_loss(x),empty)[1])
    @test_throws ArgumentError differentiable_persistence(X;modulus=3)
    @test_throws ArgumentError differentiable_persistence(X;dim=2)
    @test_throws ArgumentError differentiable_persistence(X;input=:alpha)
    @test_throws ArgumentError topological_loss(X;power=0.5)
    @test_throws ArgumentError topological_loss(X;target=-1)
    @test_throws ArgumentError differentiable_persistence([1.0 2;2 1];input=:distances)
end

# Run the same file with Enzyme loaded to verify the optional LLVM backend.
if Base.find_package("Enzyme")!==nothing
    @eval using Enzyme
    @testset "Enzyme Rips gradients" begin
        X = [0.0 1.0 1.23 0.11; 0.0 0.07 1.1 1.02]
        f = x -> topological_loss(x;target=0.7)
        gradient = zeros(size(X))
        Enzyme.autodiff(Enzyme.Reverse,f,Enzyme.Active,Enzyme.Duplicated(X,gradient))
        @test gradient ≈ central_gradient(f,X) atol=1e-6 rtol=1e-5
        for column in (1,2)
            endpoint = x -> sum(differentiable_persistence(x)[:,column])
            grad_endpoint = zeros(size(X))
            Enzyme.autodiff(Enzyme.Reverse,endpoint,Enzyme.Active,Enzyme.Duplicated(X,grad_endpoint))
            @test grad_endpoint ≈ central_gradient(endpoint,X) atol=1e-6 rtol=1e-5
        end
        D = [norm(X[:,i]-X[:,j]) for i in axes(X,2),j in axes(X,2)]
        g = d -> topological_loss(d;input=:distances)
        grad_D = zeros(size(D))
        Enzyme.autodiff(Enzyme.Reverse,g,Enzyme.Active,Enzyme.Duplicated(D,grad_D))
        @test grad_D ≈ Zygote.gradient(g,D)[1] atol=1e-6
        empty = X[:,1:2]
        grad_empty = zeros(size(empty))
        Enzyme.autodiff(Enzyme.Reverse,f,Enzyme.Active,Enzyme.Duplicated(empty,grad_empty))
        @test all(iszero,grad_empty)
    end
end
