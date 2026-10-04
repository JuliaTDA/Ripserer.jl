using SafeTestsets
using Test

@safetestset "primefield" begin
    include("base/primefield.jl")
end
@safetestset "simplex" begin
    include("base/simplex.jl")
end
@safetestset "chainelement" begin
    include("base/chainelement.jl")
end
@safetestset "chain" begin
    include("base/chain.jl")
end
@safetestset "plotting" begin
    include("base/plotting.jl")
end

@safetestset "reducedmatrix" begin
    include("computation/reducedmatrix.jl")
end
@safetestset "zerodimensional" begin
    include("computation/zerodimensional.jl")
end
@safetestset "zigzag" begin
    include("computation/zigzag.jl")
end
@safetestset "vineyards" begin
    include("computation/vineyards.jl")
end

@safetestset "rips" begin
    include("filtrations/rips.jl")
end
@safetestset "cubical" begin
    include("filtrations/cubical.jl")
end
@safetestset "custom" begin
    include("filtrations/custom.jl")
end
@safetestset "alpha" begin
    include("filtrations/alpha.jl")
end
@safetestset "edgecollapse" begin
    include("filtrations/edgecollapse.jl")
end
@safetestset "witness" begin
    include("filtrations/witness.jl")
end
@safetestset "new filtrations" begin
    include("filtrations/newfiltrations.jl")
end

@safetestset "cycles" begin
    include("extra/cycles.jl")
end
@safetestset "circularcoordinates" begin
    include("extra/circularcoordinates.jl")
end
@safetestset "mlj" begin
    include("extra/mlj.jl")
end

@safetestset "differentiable persistence forward" begin
    include("extra/differentiable_forward.jl")
end

@safetestset "P2 workflows" begin
    include("extra/p2_workflows.jl")
end

if VERSION>=v"1.9"
    @safetestset "P2 optional extensions" begin
        include("extra/p2_extensions.jl")
    end
end

@safetestset "aqua" begin
    include("aqua.jl")
end
@safetestset "doctests" begin
    include("doctests.jl")
end
