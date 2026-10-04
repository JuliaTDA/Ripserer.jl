using TDARipserer, Test, LinearAlgebra
@testset "differentiable persistence forward API" begin
    X = [0.0 1.0 1.23 0.11; 0.0 0.07 1.1 1.02]
    for dim in (0,1)
        endpoints = differentiable_persistence(X;dim=dim)
        diagram = filter(isfinite,ripserer([Tuple(X[:,i]) for i in axes(X,2)];dim_max=dim)[dim+1])
        @test endpoints[:,1] ≈ birth.(diagram)
        @test endpoints[:,2] ≈ death.(diagram)
        @test differentiable_persistence(X,persistence_pairing(X;dim=dim))==endpoints
    end
    @test size(differentiable_persistence(X[:,1:1];dim=0))==(0,2)
    @test size(differentiable_persistence(X[:,1:2];dim=1))==(0,2)
    @test topological_loss(X[:,1:2])==0
    @test_throws ArgumentError differentiable_persistence(X;modulus=3)
    @test_throws ArgumentError differentiable_persistence(X;dim=2)
end
