using Test
using TDARipserer
using TDAPersistenceDiagrams: PersistenceInterval
using JuMP, HiGHS, DelayEmbeddings

@testset "Certified optimal cycles and volumes over F₂" begin
    f=Rips([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
    interval=only(ripserer(f;reps=true)[2])
    cycle=optimal_cycle(f,interval)
    @test length(cycle.chain)==cycle.objective==4
    @test optimal_cycle(f,interval).chain==cycle.chain
    volume=optimal_volume(f,interval)
    @test length(volume.chain)==volume.objective==2
    cells=TDARipserer._chain_cells(f,2,death(interval),100)
    B=TDARipserer._mod2_boundary(cells[2],cells[3])
    c=TDARipserer._seed_vector(cells[2],cycle.chain)
    z=TDARipserer._seed_vector(cells[3],volume.chain)
    @test all(iseven,B*z-c)
    # Independent exhaustive search of all possible fillings.
    costs=[sum(v) for bits in 0:2^length(cells[3])-1
        for v in [[Int((bits>>j)&1) for j in 0:length(cells[3])-1]] if all(iseven,B*v-c)]
    @test minimum(costs)==volume.objective

    custom=Custom([(1,2)=>0.,(2,3)=>0.,(3,4)=>0.,(1,4)=>0.,
        (1,3)=>1.,(1,2,3)=>1.,(1,3,4)=>3.])
    int=only(ripserer(custom;alg=:homology,reps=true)[2])
    shorter=optimal_cycle(custom,int;at=1.,representative_kind=:cycle)
    @test shorter.objective==3
    seed=TDARipserer._interval_chain(custom,int,1.,:cycle)
    cs=TDARipserer._chain_cells(custom,2,1.,100)
    b=TDARipserer._mod2_boundary(cs[2],cs[3])
    cc=TDARipserer._seed_vector(cs[2],seed)
    zz=TDARipserer._seed_vector(cs[2],shorter.chain)
    yy=TDARipserer._seed_vector(cs[3],shorter.filling)
    @test all(iseven,zz-cc-b*yy)
    @test shorter.objective==minimum(sum(mod.(cc+b*v,2)) for bits in 0:2^size(b,2)-1
        for v in [[Int((bits>>j)&1) for j in 0:size(b,2)-1]])

    tetra=Custom([(1,2,3)=>0.,(1,2,4)=>0.,(1,3,4)=>0.,(2,3,4)=>0.,(1,2,3,4)=>1.])
    h2=only(ripserer(tetra;dim_max=2,alg=:homology,reps=true)[3])
    @test optimal_cycle(tetra,h2;representative_kind=:cycle).objective==4
    @test optimal_volume(tetra,h2;representative_kind=:cycle).objective==1
    @test_throws ArgumentError optimal_cycle(f,interval;at=death(interval))
    @test_throws ArgumentError optimal_cycle(f,interval;weights=x -> -1)
    @test_throws ArgumentError optimal_cycle(f,interval;max_cells=1)
    @test_throws ErrorException optimal_volume(f,cycle.chain;at=birth(interval))
end

@testset "DelayEmbeddings composition" begin
    s=sin.(range(0,20π;length=300))
    X=delay_embedding(s;dimension=3,delay=4)
    @test length(X)==length(s)-8
    @test collect(X)==collect(DelayEmbeddings.embed(s,3,4))
    @test first(X)≈s[[1,5,9]]
    @test suggest_delay(s;method="ac_zero",lags=1:30) in 1:30
    @test length(sliding_window_ph(X;window=20,step=30))==10
    @test_throws ArgumentError delay_embedding(s;dimension=100,delay=10)
    @test_throws ArgumentError suggest_delay(ones(20))
end
