using Test
using TDARipserer

@testset "Window PH and finite image thresholds" begin
    X=[[cos(t),sin(t)] for t in range(0,2π;length=25)[1:24]]
    windows=sliding_window_ph(X;window=12,step=6)
    @test length(windows)==3
    @test [w.start for w in windows]==[1,7,13]
    @test windows[1].stop==12
    @test windows[1].diagrams==ripserer(X[1:12])
    @test_throws ArgumentError sliding_window_ph(X;window=25)
    @test_throws ArgumentError sliding_window_ph(X;window=12,step=0)
    A=[0. 0. 0.;0. Inf 0.;0. 0. 0.]
    D=image_persistence(A)
    @test length(D[2])==1
    @test !isfinite(D[2][1])
    @test_throws ArgumentError image_persistence(fill(Inf,2,2))
    @test_throws ArgumentError image_persistence([NaN 0.;0. 0.])
end
