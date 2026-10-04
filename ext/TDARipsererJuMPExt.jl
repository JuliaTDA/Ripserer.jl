module TDARipsererJuMPExt
using TDARipserer
using SparseArrays
using JuMP
using HiGHS

function TDARipserer._solve_mod2(B::SparseMatrixCSC,c,w,mode;tie_break=:lexicographic,
        tolerance=1e-8,optimizer=HiGHS.Optimizer,time_limit=Inf)
    tie_break in (:lexicographic,:solver) || throw(ArgumentError("invalid tie_break"))
    isfinite(tolerance) && tolerance>=0 || throw(ArgumentError("invalid tolerance"))
    time_limit>0 || throw(ArgumentError("positive time_limit required"))
    model=Model(optimizer);set_silent(model)
    isfinite(time_limit) && set_time_limit_sec(model,time_limit)
    n=length(w)
    @variable(model,z[1:n],Bin)
    @variable(model,q[1:length(c)],Int)
    if mode==:homologous
        @variable(model,y[1:size(B,2)],Bin)
        @constraint(model,[i=1:length(c)],z[i]-c[i]-sum(B[i,j]*y[j] for j in 1:size(B,2))==2q[i])
    else
        y=VariableRef[]
        @constraint(model,[i=1:length(c)],sum(B[i,j]*z[j] for j in 1:n)-c[i]==2q[i])
    end
    cost=@expression(model,sum(w[j]*z[j] for j in 1:n))
    @objective(model,Min,cost)
    function solve!()
        optimize!(model)
        termination_status(model)==MOI.OPTIMAL || throw(ErrorException("chain optimization failed: $(termination_status(model)); no certified optimum"))
    end
    solve!()
    optimum=objective_value(model)
    if tie_break==:lexicographic
        @constraint(model,cost<=optimum+tolerance)
        for j in 1:n
            @objective(model,Min,z[j])
            solve!()
            fix(z[j],round(Int,value(z[j]));force=true)
        end
        solve!()
    end
    zz=Bool[round(Int,value(v))==1 for v in z]
    yy=Bool[round(Int,value(v))==1 for v in y]
    # Check the parity certificate independently of the solver's feasibility tolerance.
    residual=mode==:homologous ? Int.(zz)-c-B*Int.(yy) : B*Int.(zz)-c
    all(iseven,residual) || error("invalid integer chain certificate")
    return zz,yy,sum(w[zz])
end
end
