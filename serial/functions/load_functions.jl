# */ --------------------------------------------------------------------------------------------- /* #
# */ --- help functions -------------------------------------------------------------------------- /* #
# */ --------------------------------------------------------------------------------------------- /* #

function circ(h::Int64,H::UnitRange{Int64})::Int64
    
    (h==minimum(H)) ? idx = maximum(H) : idx = h-1
    
    return idx
end

# */ --------------------------------------------------------------------------------------------- /* #
# */ --- optimization models --------------------------------------------------------------------- /* #
# */ --------------------------------------------------------------------------------------------- /* #

function MP!(m::JuMP.Model,ms::ms_type,mp::mp_type,ps::ps_type,pp::pp_type,unc::u_type)::JuMP.Model
    
    # problem structure
    # min   f + ∑ᵢ πᵢ g(xᵢ,cᵢ)
    
    xr,cr=1:unc.nx,1:unc.nc
    # */ -- master variables --------------------------------------- /* #
    @variable(m, f, container=Array) # investment-only cost
    @variable(m, x0[ms.P,ms.I0] >= .0, container=Array) # newly installed capacity
    @variable(m, x[xr,ms.I], container=Array) # rhs paramters for operational subproblem (eg, accumulated capacity)
    @variable(m, β[ms.I] >= .0, container=Array) # operational cost of node i
    # */ -- subproblems variables ---------------------------------- /* #
    @variable(m, ϕ[ms.I,cr], container=Array) # operational cost dependent on uncertain cost cᵢ
    @variable(m, c0[ms.I], container=Array) # operational cost independent of uncertain costs c
    @variable(m, yG[ms.I,ps.G,ps.S,ps.H] >= 0, container=Array) # generation level of conventional generation
    @variable(m, yI[ms.I,ps.B,ps.S,ps.H] >= 0, container=Array) # charging power of storage generation 
    @variable(m, yO[ms.I,ps.B,ps.S,ps.H] >= 0, container=Array) # discharging power of storage generation
    @variable(m, yL[ms.I,ps.B,ps.S,ps.H] >= 0, container=Array) # energy level of storage generation
    @variable(m, yS[ms.I,ps.S,ps.H] >= 0, container=Array) # load shedding
    @variable(m, x1[ms.I,xr], container=Array) # values of rhs parameters in the subproblem (eg, capacity)
    # */ ----------------------------------------------------------- /* #
    
    # */ -- obj function ------------------------------------------- /* #
    @objective(m, Min, f + mp.κ*sum(mp.π[i]*β[i] for i in ms.I) ) # investment plus operational cost (10⁶£)
    # */ ----------------------------------------------------------- /* #
    
    # */ -- master constraints ------------------------------------- /* #
    @constraint(m, cm01, f >= exp10(-6)*sum(mp.π0[i0]*sum(mp.ci[p,i0]*x0[p,i0] for p in ms.P) for i0 in ms.I0) + exp10(-6)*mp.κ*sum(mp.π[i]*sum(mp.cf[p]*x[p,i] for p in ms.P) for i in ms.I) ) # compute investment-only cost
    @constraint(m, cm02[p=ms.P,i=ms.I], x[p,i] == mp.xh[p,i] + sum(x0[p,i0] for i0 in mp.map[i]) ) # compute accumulated capacity at node i
    @constraint(m, cm03[p=ms.P,i=ms.I], x[p,i] >= .0 ) # minimum accumulated capacity 
    @constraint(m, cm04[p=ms.P,i=ms.I], x[p,i] <= mp.xm[p] ) # maximum accumulated capacity 
    for i in ms.I for j in 1:unc.nh
        fix(x[unc.nx0+j,i],unc.h[i,j];force=true) # set rhs paramters to dependent on subproblem i (eg, CO₂ budget parameter)
    end end
    # */ -- subproblems constraints -------------------------------- /* #
    @constraint(m, cs00[i=ms.I], β[i] >= c0[i] + unc.c[i,1]*ϕ[i,1] + unc.c[i,2]*ϕ[i,2]) # compute generation cost (10⁶£) of subproblem i
    @constraint(m, cs01[i=ms.I], c0[i] == exp10(-6)*sum(sum(sum((pp.cvg[g]+pp.cfg[g]/pp.ηg[g])*yG[i,g,s,h] for g in ps.G0) + pp.cvg[ps.gn]*yG[i,ps.gn,s,h] + pp.cs*yS[i,s,h] for h in ps.H)*pp.α for s in ps.S)) # compute operational cost independent of uncertain costs c 
    @constraint(m, cs02[i=ms.I], ϕ[i,1] == exp10(-6)*sum(sum(sum(pp.eg[g]/pp.ηg[g]*yG[i,g,s,h] for g in ps.G) for h in ps.H)*pp.α for s in ps.S)) # compute cost dependent on CO₂ emission cost
    @constraint(m, cs03[i=ms.I], ϕ[i,2] == exp10(-6)*sum(sum(1/pp.ηg[ps.gn]*yG[i,ps.gn,s,h] for h in ps.H)*pp.α for s in ps.S)) # compute cost dependent on nuclear fuel cost
    @constraint(m, cs04[i=ms.I,b=ps.B,s=ps.S,h=ps.H], yL[i,b,s,h] - yL[i,b,s,circ(h,ps.H)] == pp.ηb[b]*yI[i,b,s,h] - yO[i,b,s,h]) # storage energy balance
    @constraint(m, cs05[i=ms.I,g=ps.G,s=ps.S,h=ps.H], yG[i,g,s,h] - yG[i,g,s,circ(h,ps.H)] <= pp.rg[g]*x1[i,g]) # upward ramping limitation on conventional generation
    @constraint(m, cs06[i=ms.I,g=ps.G,s=ps.S,h=ps.H], yG[i,g,s,circ(h,ps.H)] - yG[i,g,s,h] <= pp.rg[g]*x1[i,g]) # downward ramping limitation on conventional generation
    @constraint(m, cs07[i=ms.I,g=ps.G,s=ps.S,h=ps.H], yG[i,g,s,h] <= x1[i,g]) # maximum capacity limitation on conventional generation
    @constraint(m, cs08[i=ms.I,b=ps.B,s=ps.S,h=ps.H], yI[i,b,s,h] <= pp.pb[b]*x1[i,b+ps.b0]) # maximum capacity limitation on storage charging power
    @constraint(m, cs09[i=ms.I,b=ps.B,s=ps.S,h=ps.H], yO[i,b,s,h] <= pp.pb[b]*x1[i,b+ps.b0]) # maximum capacity limitation on storage discharging power
    @constraint(m, cs10[i=ms.I,b=ps.B,s=ps.S,h=ps.H], yL[i,b,s,h] <= x1[i,b+ps.b0]) # maximum capacity limitation on storage energy level
    @constraint(m, cs11[i=ms.I], sum(yG[i,g,s,h]*pp.eg[g]/pp.ηg[g] for g in ps.G for s in ps.S for h in ps.H) <= pp.lco*x1[i,unc.nx0+1]) # CO₂ emission budget
    @constraint(m, cs12[i=ms.I,s=ps.S,h=ps.H], sum(yG[i,g,s,h] for g in ps.G) + sum(yO[i,b,s,h]-yI[i,b,s,h] for b in ps.B) + yS[i,s,h] >= -x1[i,unc.nx0+2]*pp.pd[s,h] - sum(pp.pr[r,s,h]*x1[i,r+ps.r0] for r in ps.R)) # energy demand
    @constraint(m, cs13[i=ms.I,n=1:unc.nx0], x1[i,n]*exp10(-3) == x[n,i]) # generation to level of decisions {xᵢ, ∀i} set by the master problem
    @constraint(m, cs14[i=ms.I,n=1:unc.nh ], x1[i,unc.nx0+n]   == x[unc.nx0+n,i]) # rhs uncertain parameters (energy level and CO₂ budget) specific to problem i
    # */ ----------------------------------------------------------- /* #
    
    return m
end

function RMP!(m::JuMP.Model,ms::ms_type,mp::mp_type,unc::u_type)::JuMP.Model
    
    # problem structure
    # min   f + ∑ᵢ πᵢ βᵢ)
    # s.t.  f = q₀ᵀx₀ + ∑ᵢ qᵢᵀxᵢ
    #      (x₀,x₁,..,xᵢ,..,xₙ) ∈ 𝕏
    #       βᵢ ∈ Θᵢ,  ∀i
    
    xr=1:unc.nx
    # */ -- variables ---------------------------------------------- /* #
    @variable(m, f, container=Array) # investment-only cost
    @variable(m, x0[ms.P,ms.I0] >= .0, container=Array) # newly installed capacity
    @variable(m, x[xr,ms.I], container=Array) # rhs paramters for operational subproblem (eg, accumulated capacity)
    @variable(m, β[ms.I] >= .0, container=Array) # operational cost of node i
    # */ ----------------------------------------------------------- /* #
    
    # */ -- obj function ------------------------------------------- /* #
    @objective(m, Min, f + mp.κ*sum(mp.π[i]*β[i] for i in ms.I) ) # investment plus operational cost (10⁶£)
    # */ ----------------------------------------------------------- /* #
    
    # */ -- constraints -------------------------------------------- /* #
    @constraint(m, c01, f >= exp10(-6)*sum(mp.π0[i0]*sum(mp.ci[p,i0]*x0[p,i0] for p in ms.P) for i0 in ms.I0) + exp10(-6)*mp.κ*sum(mp.π[i]*sum(mp.cf[p]*x[p,i] for p in ms.P) for i in ms.I) ) # compute investment-only cost
    @constraint(m, c02[p=ms.P,i=ms.I], x[p,i] == mp.xh[p,i] + sum(x0[p,i0] for i0 in mp.map[i]) ) # compute accumulated capacity at node i
    @constraint(m, c03[p=ms.P,i=ms.I], x[p,i] >= .0 ) # minimum accumulated capacity 
    @constraint(m, c04[p=ms.P,i=ms.I], x[p,i] <= mp.xm[p] ) # maximum accumulated capacity 
    for i in ms.I for j in 1:unc.nh
        fix(x[unc.nx0+j,i],unc.h[i,j];force=true) # set rhs paramters to dependent on subproblem i (eg, CO₂ budget parameter)
    end end
    # */ ----------------------------------------------------------- /* #
    optimize!(m)
    
    return m
end

function SP!(m::JuMP.Model,ps::ps_type,pp::pp_type,unc::u_type)::JuMP.Model
    
    # problem structure
    # min   c₀ + cᵀϕ
    # s.t.  c₀  = C₀y
    #       ϕ   = C y
    #       y ∈ 𝕐
    #       A y ≦ B x
    #       x = xf :(λ)
    
    # */ -- variables ---------------------------------------------- /* #
    @variable(m, ϕ[1:unc.nc], container=Array) # operational cost dependent on uncertain cost cᵢ
    @variable(m, c0, container=Array) # operational cost independent of uncertain costs c
    @variable(m, yG[ps.G,ps.S,ps.H] >= 0, container=Array) # generation level of conventional generation
    @variable(m, yI[ps.B,ps.S,ps.H] >= 0, container=Array) # charging power of storage generation
    @variable(m, yO[ps.B,ps.S,ps.H] >= 0, container=Array) # discharging power of storage generation
    @variable(m, yL[ps.B,ps.S,ps.H] >= 0, container=Array) # energy level of storage generation
    @variable(m, yS[ps.S,ps.H] >= 0, container=Array) # load shedding
    @variable(m, x0[1:unc.nx], container=Array) # values of rhs parameters in the subproblem (eg, capacity)
    @variable(m, x[1:unc.nx], container=Array) # values of decisions {xᵢ, ∀i} set by the master problem
    # */ ----------------------------------------------------------- /* #
    
    # */ -- obj function ------------------------------------------- /* #
    @objective(m, Min, c0 ) # generation cost (10⁶£)
    # */ ----------------------------------------------------------- /* #
    
    # */ -- constraints -------------------------------------------- /* #
    @constraint(m, c01, c0 == exp10(-6)*sum(sum(sum((pp.cvg[g]+pp.cfg[g]/pp.ηg[g])*yG[g,s,h] for g in ps.G0) + pp.cvg[ps.gn]*yG[ps.gn,s,h] + pp.cs*yS[s,h] for h in ps.H)*pp.α for s in ps.S)) # compute operational cost independent of uncertain costs c 
    @constraint(m, c02, ϕ[1] == exp10(-6)*sum(sum(sum(pp.eg[g]/pp.ηg[g]*yG[g,s,h] for g in ps.G) for h in ps.H)*pp.α for s in ps.S)) # compute cost dependent on CO₂ emission cost
    @constraint(m, c03, ϕ[2] == exp10(-6)*sum(sum(1/pp.ηg[ps.gn]*yG[ps.gn,s,h] for h in ps.H)*pp.α for s in ps.S)) # compute cost dependent on nuclear fuel cost
    @constraint(m, c04[b=ps.B,s=ps.S,h=ps.H], yL[b,s,h] - yL[b,s,circ(h,ps.H)] == pp.ηb[b]*yI[b,s,h] - yO[b,s,h]) # storage energy balance
    @constraint(m, c05[g=ps.G,s=ps.S,h=ps.H], yG[g,s,h] - yG[g,s,circ(h,ps.H)] <= pp.rg[g]*x0[g]) # upward ramping limitation on conventional generation
    @constraint(m, c06[g=ps.G,s=ps.S,h=ps.H], yG[g,s,circ(h,ps.H)] - yG[g,s,h] <= pp.rg[g]*x0[g]) # downward ramping limitation on conventional generation
    @constraint(m, c07[g=ps.G,s=ps.S,h=ps.H], yG[g,s,h] <= x0[g]) # maximum capacity limitation on conventional generation
    @constraint(m, c08[b=ps.B,s=ps.S,h=ps.H], yI[b,s,h] <= pp.pb[b]*x0[b+ps.b0]) # maximum capacity limitation on storage charging power
    @constraint(m, c09[b=ps.B,s=ps.S,h=ps.H], yO[b,s,h] <= pp.pb[b]*x0[b+ps.b0]) # maximum capacity limitation on storage discharging power
    @constraint(m, c10[b=ps.B,s=ps.S,h=ps.H], yL[b,s,h] <= x0[b+ps.b0]) # maximum capacity limitation on storage energy level
    @constraint(m, c11, sum(yG[g,s,h]*pp.eg[g]/pp.ηg[g] for g in ps.G for s in ps.S for h in ps.H) <= pp.lco*x0[unc.nx0+1]) # CO₂ emission budget
    @constraint(m, c12[s=ps.S,h=ps.H], sum(yG[g,s,h] for g in ps.G) + sum(yO[b,s,h]-yI[b,s,h] for b in ps.B) + yS[s,h] >= -x0[unc.nx0+2]*pp.pd[s,h] - sum(pp.pr[r,s,h]*x0[r+ps.r0] for r in ps.R)) # energy demand
    @constraint(m, c13[n=1:unc.nx0], x0[n]*exp10(-3) == x[n]) # generation to level of decisions {xᵢ, ∀i} set by the master problem
    @constraint(m, c14[n=1:unc.nh ], x0[unc.nx0+n]   == x[unc.nx0+n]) # rhs uncertain parameters (energy level and CO₂ budget) specific to problem i
    # */ ----------------------------------------------------------- /* #
    optimize!(m)
    
    return m
end

function LBoracle!(m::JuMP.Model,unc::u_type,J::Int64)::JuMP.Model
    
    # problem structure
    # min   ϕ + γᵀcᵢ
    # s.t.  ϕ + γᵀcₛ ≧ θₛ + λᵀₛ(x-xₛ),   ∀s
    #       x = xᵢ
    
    # */ -- variables ---------------------------------------------- /* #
    @variable(m, ϕ, container=Array)
    @variable(m, γ[1:unc.nc] >= .0, container=Array)
    @variable(m, x[1:unc.nx], container=Array)
    # */ ----------------------------------------------------------- /* #
    
    # */ -- obj function ------------------------------------------- /* #
    @objective(m, Min, 0. )
    # */ ----------------------------------------------------------- /* #
    
    # */ -- constraints -------------------------------------------- /* #
    # */ ----------------------------------------------------------- /* #
    optimize!(m)
    
    return m
end

function UBoracle!(m::JuMP.Model,unc::u_type,J::Int64)::JuMP.Model
    
    # problem structure
    # max   ϕ + γᵀxᵢ
    # s.t.  ϕ + γᵀxₛ ≦ θₛ + ϕᵀₛ(c-cₛ),   ∀s
    #       c = cᵢ
    
    # */ -- variables ---------------------------------------------- /* #
    @variable(m, ϕ, container=Array)
    @variable(m, γ[1:unc.nx] <= .0, container=Array)
    @variable(m, c[1:unc.nc], container=Array)
    # */ ----------------------------------------------------------- /* #
    
    # */ -- obj function ------------------------------------------- /* #
    @objective(m, Max, 0. )
    # */ ----------------------------------------------------------- /* #
    
    # */ -- constraints -------------------------------------------- /* #
    # */ ----------------------------------------------------------- /* #
    optimize!(m)
    
    return m
end

# */ --------------------------------------------------------------------------------------------- /* #
# */ --- functions to generate structs : master -------------------------------------------------- /* #
# */ --------------------------------------------------------------------------------------------- /* #

function gen_kiidx(u::u_type,w::Int64)::Tuple{Array{BitArray{1},1},Array{Array{Int64,1},1}}

    # function to generate (kidx,iidx) for given number of subproblems w solved at each iteration
    # iidx -> array of subsets {ℰₖ, 𝑘=1,..,w}
    # kidx -> array of 0/1 association with subset ℰₖ 

    if w == 1 # (trivial case with one subproblem) ℰₖ := ℰ
        kidx = Array{BitArray{1},1}(undef,0)
        push!(kidx,zeros(Int64,u.ni).==0)
        iidx = Array{Array{Int64,1},1}(undef,0)
        push!(iidx,[i for i in 1:u.ni])
    else
        x = zeros(u.nh+u.nc,u.ni) # generate empty matrix x for (H,C) values of each subproblem
        for j in 1:u.nh
            uh = u.h[:,j] # collect values of H (rhs uncertain params)
            hmin,hmax = minimum(uh),maximum(uh) # evaluate min & max values 
            (hmax > hmin) ? x[j,:] .= (uh.-hmin)./(hmax-hmin) :  x[j,:] .= zeros(u.ni) # set normalized H values into x
        end
        for j in 1:u.nc
            uc = u.c[:,j] # collect values of C (uncertain cost coeff)
            cmin,cmax = minimum(uc),maximum(uc) # evaluate min & max values 
            (cmax > cmin) ? x[j+u.nh,:] .= (uc.-cmin)./(cmax-cmin) : x[j+u.nh,:] .= zeros(u.ni) # set normalized C values into x
        end
        kass = kmeans(x,w).assignments # run kmeans clustering on x to obtain w subsets
        kidx = Array{BitArray{1},1}(undef,0) # generate empty kidx
        for j in 1:w 
            push!(kidx,kass.==j) # for each subset j = 1,..,w generate a bitarry with 1 if subproblem i belongs to subset j (and 0 otherwise)
        end
        idx = [i for i in 1:u.ni] # generate array of indices of subproblems
        iidx = Array{Array{Int64,1},1}(undef,0) # generate empty iidx
        for j in 1:w
            push!(iidx,idx[kidx[j]]) # for each subset j = 1,..,w generate an array with all the subproblems that belongs to subset j
        end
    end

    return kidx,iidx
end

function gen_D__type(ms::ms_type,mp::mp_type,ps::ps_type,pp::pp_type,u::u_type,c::Int64,a::Int64,q::Float64,w::Int64,j::Int64,e::Float64)::D__type

    # function to generate D__type (structure of problem data) for given ms (master problem sets), mp (master problem params), 
    # ps (subproblem sets), pp (subproblem params), u (uncertainty params), c (case), a (algorithm), q (tolerance damping),
    # w (number of exact subproblems), j (maximum number of iteration), and e (convergence tolerance)
    
    k,i = gen_kiidx(u,w) # generate (kidx,iidx) for given uncertain data u and given number of subproblems w solved at each iteration

    return D__type(ms,mp,ps,pp,u,c,a,q,w,k,i,j,e)
end

function gen_R__type(d::D__type)::R__type

    # function to generate R__type (structure of relax master problem) given data D__type (problem data)
    
    m = JuMP.Model(JuMP.optimizer_with_attributes(()->Gurobi.Optimizer(gurobi_env))) # generate empty JuMP optimization model
    JuMP.set_optimizer_attribute(m,"OutputFlag",0) # set JuMP model silent
    RMP!(m,d.ms,d.mp,d.unc) # build relaxed master problem
    v = Array{JuMP.VariableRef,2}(undef,d.unc.nx+1,d.unc.ni) # generate empty array of vars
    for i in 1:d.unc.ni v[:,i] .= vcat([m[:β][i]],m[:x][:,i]) end # set reference to variablems {(βᵢ,xᵢ), ∀i}
    c = Array{JuMP.ConstraintRef,2}(undef,d.J,d.unc.ni) # generate empty array of cons
    
    return R__type(m,v,c)
end

function gen_T1_type(d::D__type)::T1_type

    # function to generate T1_type (structure of temporary data, algorithm 1 ) given data D__type (problem data)
    
    x = zeros(d.unc.nx,d.unc.ni) # generate empty array for decisions {xᵢ, ∀i}
    c = convert(Array{Float64,2},d.unc.c') # set values of parameters {cᵢ, ∀i} 
    β = zeros(d.unc.ni) # generate empty array for variables value {βᵢ, ∀i}
    θ = zeros(d.unc.ni) # generate empty array for optimal objectives {θᵢ, ∀i}
    λ = zeros(d.unc.nx,d.unc.ni) # generate empty array for subgradient {λᵢ, ∀i} wrt xᵢ 
    
    return T1_type(x,c,β,θ,λ,0.)
end

function gen_T3_type(d::D__type)::T3_type

    # function to generate T1_type (structure of temporary data, algorithm 3) given data D__type (problem data)
    
    x  = zeros(d.unc.nx,d.unc.ni) # generate empty array for decisions {xᵢ, ∀i}
    c  = convert(Array{Float64,2},d.unc.c') # set values of parameters {cᵢ, ∀i}
    β  = zeros(d.unc.ni) # generate empty array for variables value {βᵢ, ∀i}
    θp = zeros(d.unc.ni) # generate empty array for valid upper bounds on objectives {θᵢ, ∀i}
    θd = zeros(d.unc.ni) # generate empty array for valid lower bounds on objectives {θᵢ, ∀i}
    λ  = zeros(d.unc.nx,d.unc.ni) # generate empty array for valid subgradient {λᵢ, ∀i} wrt xᵢ 
     
    return T3_type(x,c,β,θp,θd,λ,0.)
end

function gen_T2_type(d::D__type)::T2_type

    # function to generate T1_type (structure of temporary data, algorithm 2) given data D__type (problem data)
    
    x  = zeros(d.unc.nx,d.unc.ni) # generate empty array for decisions {xᵢ, ∀i}
    c  = convert(Array{Float64,2},d.unc.c') # set values of parameters {cᵢ, ∀i}
    β  = zeros(d.unc.ni) # generate empty array for variables value {βᵢ, ∀i}
    θl = zeros(d.unc.ni) # generate empty array for valid lower bounds on objectives {θᵢ, ∀i}
    λl = zeros(d.unc.nx,d.unc.ni) # generate empty array for valid subgradient {λᵢ, ∀i} wrt xᵢ 
    θu = zeros(d.unc.ni) # generate empty array for valid upper bounds on objectives {θᵢ, ∀i}
    ϕu = zeros(d.unc.nc,d.unc.ni) # generate empty array for valid subgradient {ϕᵢ, ∀i} wrt cᵢ
    Ie = zeros(Int64,d.w) # generate empty array for indices of w subproblems
    
    return T2_type(x,c,β,θl,λl,θu,ϕu,Ie,0.)
end

function gen_H__type(d::D__type)::H__type

    # function to generate H_type (structure of stored data) given data D__type (problem data)
    
    L = zeros(d.J) # generate empty array for lower bound values
    U = zeros(d.J) # generate empty array for upper bound values
    t = zeros(d.J) # generate empty array for time spend to perform each iteration
    T = zeros(d.J,6) # generate empty array for time spend to perform each iteration, split by steps
    
    return H__type(L,U,t,T,0)
end

function gen_M__type(d::D__type)::M__type

    # function to generate M_type (structure of exact solutions stored) given data D__type (problem data)

    θ = zeros(nprocs()) # generate empty array of optimal objective
    λ = zeros(nprocs(),d.unc.nx) # generate empty array of subgradient wrt x
    ϕ = zeros(nprocs(),d.unc.nc) # generate empty array of subgradient wrt c
    x = zeros(nprocs(),d.unc.nx) # generate empty array of values of x
    c = zeros(nprocs(),d.unc.nc) # generate empty array of values of c

    return M__type(0,θ,λ,ϕ,x,c)
end

function gen_B0_type(c::Int64)::JuMP.Model

    # function to generate B0_type (Benders master structure of algorithm 0) given c (case)
    
    ms,mp,ps,pp,unc = load_data(c) # generate data structures ms,mp,ps,pp,unc for given value of c (case)
    m = JuMP.Model(JuMP.optimizer_with_attributes(()->Gurobi.Optimizer(gurobi_env))) # generate empty JuMP model
    println(" Attempting to build JuMP model of deterministic equivalent")
    println(" ")
    MP!(m,ms,mp,ps,pp,unc) # build master problem
    println(" JuMP Model build successfully!")
    println(" ")
    
    return m
end

function gen_B1_type(c::Int64,a::Int64,j::Int64,e::Float64)::B1_type

    # function to generate B1_type (Benders master structure of algorithm 1 ) given c (case), a (algorithm), 
    # j (maximum number of iteration), and e (convergence tolerance).
    
    ms,mp,ps,pp,unc = load_data(c) # generate data structures ms,mp,ps,pp,unc for given value of c (case)
    d = gen_D__type(ms,mp,ps,pp,unc,c,a,.0,1,j,e) # generate D__type (structure of problem data)
    r = gen_R__type(d) # generate R__type (structure of relax master problem)
    t = gen_T1_type(d) # generate T1_type (structure of temporary data)
    h = gen_H__type(d) # generate H__type (structure of stored data)
    
    return B1_type(r,t,h,d) 
end

function gen_B3_type(c::Int64,a::Int64,q::Float64,j::Int64,e::Float64)::B3_type

    # function to generate B3_type (Benders master structure of algorithm 3) given c (case), a (algorithm), q (tolerance damping),
    # j (maximum number of iteration), and e (convergence tolerance).
    
    ms,mp,ps,pp,unc = load_data(c) # generate data structures ms,mp,ps,pp,unc for given value of c (case)
    d = gen_D__type(ms,mp,ps,pp,unc,c,a,q,1,j,e) # generate D__type (structure of problem data)
    r = gen_R__type(d) # generate R__type (structure of relax master problem)
    t = gen_T3_type(d) # generate T3_type (structure of temporary data)
    h = gen_H__type(d) # generate H__type (structure of stored data)
    
    return B3_type(r,t,h,d) 
end

function gen_B2_type(c::Int64,a::Int64,w::Int64,j::Int64,e::Float64)::B2_type

    # function to generate B2_type (Benders master structure of algorithm 2) given c (case), a (algorithm), w (number of exact subproblems),
    # j (maximum number of iteration), and e (convergence tolerance).
    
    ms,mp,ps,pp,unc = load_data(c) # generate data structures ms,mp,ps,pp,unc for given value of c (case)
    d = gen_D__type(ms,mp,ps,pp,unc,c,a,.0,w,j,e) # generate D__type (structure of problem data)
    r = gen_R__type(d) # generate R__type (structure of relax master problem)
    t = gen_T2_type(d) # generate T2_type (structure of temporary data)
    h = gen_H__type(d) # generate H__type (structure of stored data)
    m = gen_M__type(d) # generate M__type (structure of exact solutions stored)
    
    return B2_type(r,t,h,m,d)
end

# */ --------------------------------------------------------------------------------------------- /* #
# */ --- functions to generate structs : subproblem ---------------------------------------------- /* #
# */ --------------------------------------------------------------------------------------------- /* #

function gen_d__type(d::D__type)::d__type

    # function to generate d__type (structure of subproblem data) for given D__type (structure of problem data) 

    s = d.ps # set subproblem sets
    p = d.pp # set subproblem parameters
    u = d.unc # set uncertainty parameters
    a = d.algm # set algorithm
    j = Int64(round(d.J*d.w;digits=0)) # set maximum number of exact solutions to store
    
    return d__type(s,p,u,a,j)
end

function gen_Lb_type(d::d__type)::O__type

    # function to generate O__type (structure of oracle for lower bound) for given d__type (structure of subproblem data)

    m = JuMP.Model(JuMP.optimizer_with_attributes(()->Gurobi.Optimizer(gurobi_env))) # generate empty JuMP model
    JuMP.set_optimizer_attribute(m,"OutputFlag",0) # set model silent
    LBoracle!(m,d.unc,d.J+1) # build lower bound oracle
    v = vcat(m[:ϕ],m[:γ]) # set reference to variables (ϕ,γ) 
    c = Array{JuMP.ConstraintRef,1}(undef,d.J+1) # generate empty reference to new constraints  
    o = v'*ones(length(v)) # generate reference for objective function
    e = v'*ones(length(v)) # generate reference for objective function
    h = zeros(d.J+1) # generate empty array for rhs values of old constraints (old exact solutions)

    return O__type(m,c,v,o,e,h,0)
end

function gen_Ub_type(d::d__type)::O__type

    # function to generate O__type (structure of oracle for upper bound) for given d__type (structure of subproblem data)

    m = JuMP.Model(JuMP.optimizer_with_attributes(()->Gurobi.Optimizer(gurobi_env))) # generate empty JuMP model
    JuMP.set_optimizer_attribute(m,"OutputFlag",0) # set model silent
    UBoracle!(m,d.unc,d.J+1) # build upper bound oracle
    v = vcat(m[:ϕ],m[:γ]) # set reference to variables (ϕ,γ)
    c = Array{JuMP.ConstraintRef,1}(undef,d.J+1) # generate empty reference to new constraints  
    o = v'*ones(length(v)) # generate reference for objective function
    e = v'*ones(length(v)) # generate reference for objective function
    h = zeros(d.J+1) # generate empty array for rhs values of old constraints (old exact solutions)

    return O__type(m,c,v,o,e,h,0)
end

function gen_E__type(d::d__type)::E__type

    # function to generate E__type (structure of exact subproblem) for given d__type (structure of subproblem data)

    m = JuMP.Model(JuMP.optimizer_with_attributes(()->Gurobi.Optimizer(gurobi_env))) # generate empty JuMP model
    JuMP.set_optimizer_attribute(m,"OutputFlag",0) # set model silent
    JuMP.set_optimizer_attribute(m,"Method",2) # set solution method (barrier)
    d.algm == 3 ? JuMP.set_optimizer_attribute(m,"Crossover",0) : nothing # deactivate crossover (if algorithm 3)
    SP!(m,d.ps,d.pp,d.unc) # build subproblem
    v = vcat(m[:c0],m[:ϕ][:]) # set reference to variables (c₀,ϕ)
    o = v'*ones(length(v)) # generate reference for objective function

    return E__type(m,v,o)
end

function gen_M__type(d::d__type)::M__type

    # function to generate M__type (structure of exact solutions stored) for given d__type (structure of subproblem data)

    θ = zeros(d.J+1) # generate empty array for optimal objective of exact solutions
    λ = zeros(d.J+1,d.unc.nx) # generate empty array for subgradient wrt x of exact solutions
    ϕ = zeros(d.J+1,d.unc.nc) # generate empty array for subgradient wrt c of exact solutions
    x = zeros(d.J+1,d.unc.nx) # generate empty array for values of x of exact solutions
    c = zeros(d.J+1,d.unc.nc) # generate empty array for values of c of exact solutions

    return M__type(0,θ,λ,ϕ,x,c)
end

function gen_Q__type(d::d__type)::Array{Q__type,1}

    # function to generate Q__type (structure of extended exact solutions stored) for given d__type (structure of subproblem data)

    return Array{Q__type,1}(undef,0)
end

function gen_t1_type(d::d__type)::t1_type

    # function to generate t1_type (structure of temporary data, algorithm 1 ) for given d__type (structure of subproblem data)

    x = zeros(d.unc.nx,d.unc.ni) # generate empty array for decisions {xᵢ, ∀i}
    c = convert(Array{Float64,2},d.unc.c') # set values of parameters {cᵢ, ∀i}
    θ = zeros(d.unc.ni) # generate empty array for optimal objectives {θᵢ, ∀i}
    λ = zeros(d.unc.nx,d.unc.ni) # generate empty array for subgradient {λᵢ, ∀i} wrt xᵢ 

   return t1_type(x,c,θ,λ,0)
end

function gen_t3_type(d::d__type)::t3_type

    # function to generate t3_type (structure of temporary data, algorithm 3) for given d__type (structure of subproblem data)

    x  = zeros(d.unc.nx,d.unc.ni) # generate empty array for decisions {xᵢ, ∀i}
    c  = convert(Array{Float64,2},d.unc.c') # set values of parameters {cᵢ, ∀i}
    θp = zeros(d.unc.ni) # generate empty array for valid upper bounds on objectives {θᵢ, ∀i}
    θd = zeros(d.unc.ni) # generate empty array for valid power bounds on objectives {θᵢ, ∀i}
    λ  = zeros(d.unc.nx,d.unc.ni) # generate empty array for valid subgradients {λᵢ, ∀i} wrt xᵢ

   return t3_type(x,c,θp,θd,λ,0)
end

function gen_t2_type(d::d__type)::t2_type

    # function to generate t2_type (structure of temporary data, algorithm 2) for given d__type (structure of subproblem data)

    x  = zeros(d.unc.nx,d.unc.ni) # generate empty array for decisions {xᵢ, ∀i}
    c  = convert(Array{Float64,2},d.unc.c') # set values of parameters {cᵢ, ∀i}
    θl = zeros(d.unc.ni) # generate empty array for valid power bounds on objectives {θᵢ, ∀i}
    θu = zeros(d.unc.ni) # generate empty array for valid upper bounds on objectives {θᵢ, ∀i}
    λl = zeros(d.unc.nx,d.unc.ni) # generate empty array for valid subgradients {λᵢ, ∀i} wrt xᵢ 
    ϕu = zeros(d.unc.nc,d.unc.ni) # generate empty array for valid subgradients {ϕᵢ, ∀i} wrt cᵢ 

   return t2_type(x,c,θl,θu,λl,ϕu,0)
end

function gen_S1_type(D::D__type)::S1_type

    # function to generate S1_type (structure of Benders subproblem, algorithm 1 ) for given d__type (structure of subproblem data)

    d = gen_d__type(D) # generate d__type (structure of subproblem data) 
    e = gen_E__type(d) # generate E__type (structure of exact subproblem) 
    t = gen_t1_type(d) # generate t1_type (structure of temporary data) 

    return S1_type(e,d,t)
end

function gen_S3_type(D::D__type)::S3_type

    # function to generate S3_type (structure of Benders subproblem, algorithm 3) for given d__type (structure of subproblem data)

    d = gen_d__type(D) # generate d__type (structure of subproblem data)
    e = gen_E__type(d) # generate E__type (structure of exact subproblem) 
    t = gen_t3_type(d) # generate t3_type (structure of temporary data) 

    return S3_type(e,d,t)
end

function gen_S2_type(D::D__type)::S2_type

    # function to generate S2_type (structure of Benders subproblem, algorithm 2) for given d__type (structure of subproblem data)

    d = gen_d__type(D) # generate d__type (structure of subproblem data)
    e = gen_E__type(d) # generate E__type (structure of exact subproblem) 
    l = gen_Lb_type(d) # generate O__type (structure of oracle for lower bound)
    u = gen_Ub_type(d) # generate O__type (structure of oracle for upper bound)
    m = gen_M__type(d) # generate M__type (structure of exact solutions stored)
    q = gen_Q__type(d) # generate Q__type (structure of extended exact solutions stored)
    t = gen_t2_type(d) # generate t2_type (structure of temporary data) 

    return S2_type(e,l,u,m,q,d,t)
end

# */ --------------------------------------------------------------------------------------------- /* #
# */ --- functions to generate structs : master and subproblem ----------------------------------- /* #
# */ --------------------------------------------------------------------------------------------- /* #

function generate_structures(c::Int64,a::Int64,q::Float64,w::Int64,j::Int64,e::Float64)

    # function to generate B__struct and S__struct for given values of c (case), a (algorithm), q (tolerance damping),
    # w (number of exact subproblems), j (maximum number of iteration), and e (convergence tolerance)

    if a == 0 
        B = gen_B0_type(c) # generate B__struct (Benders master structure)
        S = nothing # empty subproblem structure
    end
    if a == 1 # (if Stand_Bend)
        B = gen_B1_type(c,a,j,e) # generate B__struct (Benders master structure)
        S = gen_S1_type(B.data) # generate S__struct (Benders subproblem structure)
    end
    if a == 2 # (if Adapt_Bend)
        B = gen_B2_type(c,a,w,j,e) # generate B__struct (Benders master structure)
        S = gen_S2_type(B.data) # generate S__struct (Benders subproblem structure)
    end
    if a == 3 # (if Zaker_Bend)
        B = gen_B3_type(c,a,q,j,e) # generate B__struct (Benders master structure)
        S = gen_S3_type(B.data) # generate S__struct (Benders subproblem structure)
    end
    
    return B,S
end

# */ --------------------------------------------------------------------------------------------- /* #
# */ --- lower level functions : subproblem ------------------------------------------------------ /* #
# */ --------------------------------------------------------------------------------------------- /* #

function set_Iex!(b::B2_type)::B2_type

    # function to update B2_type.temp.Ie (index of w subproblems solved exactly @ iter k) 
    # for given B2_type (Benders master problem) -> algorithm 2 (Adapt_Bend)

    for k in 1:b.data.w # for each k = 1,..,w (for each subset ℰₖ)
        ie = findmax(b.data.mp.π[b.data.kidx[k]].*(b.temp.θu[b.data.kidx[k]].-b.temp.θl[b.data.kidx[k]]))[2] # select subproblem ie ∈ ℰₖ for which the gap is the largest
        b.temp.Ie[k] = b.data.iidx[k][ie] # store ie in B2_type.temp.Ie
    end

    return b
end

function solv_exact!(s::S1_type)::S1_type

    # function to update and solve S1_type.ex.m (exact subproblem) and store exact solution in S1_type.temp
    # for given S1_type (Benders subproblem) -> algorithm 1 (Stand_Bend) 

    s.ex.objf =  s.ex.vars'*vcat([1.],s.temp.c[:,s.temp.i]) # compute new objective of subproblem i (c₀ + cᵢᵀϕ)
    @objective(s.ex.m, Min, s.ex.objf) # update objective in JuMP model (c₀ + cᵢᵀϕ)
    fix.(s.ex.m[:x],s.temp.x[:,s.temp.i];force=true) # fix values of variables xᵢ in the model to decisions xᵢ set by the master problem
    optimize!(s.ex.m) # solve the JuMP model to optimality
    s.temp.θ[s.temp.i] = objective_value(s.ex.m) # store optimal objective θᵢ
    s.temp.λ[:,s.temp.i] .= dual.(FixRef.(s.ex.m[:x])) # store subgradient λᵢ wrt xᵢ

    return s
end

function solv_exact!(s::S3_type)::S3_type

    # function to update and solve S3_type.ex.m (exact subproblem) and store exact solution in S3_type.temp
    # for given S3_type (Benders subproblem) -> algorithm 3 (Zaker_Bend)

    s.ex.objf =  s.ex.vars'*vcat([1.],s.temp.c[:,s.temp.i]) # compute new objective of subproblem i (c₀ + cᵢᵀϕ)
    @objective(s.ex.m, Min, s.ex.objf) # update objective in JuMP model (c₀ + cᵢᵀϕ)
    fix.(s.ex.m[:x],s.temp.x[:,s.temp.i];force=true) # fix values of variables xᵢ in the model to decisions xᵢ set by the master problem
    optimize!(s.ex.m) # solve the JuMP model up to optimality tolerance δ
    s.temp.θp[s.temp.i] = objective_value(s.ex.m) # save suboptimal primal objective θᵢ (upper bound)
    s.temp.λ[:,s.temp.i] .= dual.(FixRef.(s.ex.m[:x])) # store subgradient λᵢ wrt xᵢ
    s.temp.θd[s.temp.i] = s.temp.λ[:,s.temp.i]'*value.(s.ex.m[:x]) # save suboptimal dual objective θᵢ (lower bound)

    return s
end

function solv_exact!(s::S2_type)::S2_type

    # function to update and solve S2_type.ex.m (exact subproblem) and store exact solution in S2_type.temp
    # for given S2_type (Benders subproblem) -> algorithm 2 (Adapt_Bend)

    s.ex.objf =  s.ex.vars'*vcat([1.],s.temp.c[:,s.temp.i]) # compute new objective of subproblem i (c₀ + cᵢᵀϕ)
    @objective(s.ex.m, Min, s.ex.objf) # update objective in JuMP model (c₀ + cᵢᵀϕ)
    fix.(s.ex.m[:x],s.temp.x[:,s.temp.i];force=true) # fix values of variables xᵢ in the model to decisions xᵢ set by the master problem
    optimize!(s.ex.m) # solve the JuMP model to optimality
    s.temp.θl[s.temp.i] = objective_value(s.ex.m) # store optimal objective θᵢ
    s.temp.λl[:,s.temp.i] .= dual.(FixRef.(s.ex.m[:x])) # store subgradient λᵢ wrt xᵢ
    s.temp.ϕu[:,s.temp.i] .= value.(s.ex.m[:ϕ]) # store subgradient ϕᵢ wrt cᵢ

    return s
end

function extended_sol!(s::S2_type)::S2_type

    # function to store new extended exact solution to S2_type.q (extended exact solutions stored)
    # for given S2_type (Benders subproblem) -> algorithm 2 (Adapt_Bend)

    ϕ  = value.(s.ex.m[:ϕ] ) # extract the operational cost dependent on uncertain cost cᵢ
    c0 = value.(s.ex.m[:c0]) # extract the operational cost independent of uncertain costs c
    yG = value.(s.ex.m[:yG]) # extract the generation level of conventional generation
    yI = value.(s.ex.m[:yI]) # extract the charging power of storage generation
    yO = value.(s.ex.m[:yI]) # extract the discharging power of storage generation
    yL = value.(s.ex.m[:yL]) # extract the energy level of storage generation
    yS = value.(s.ex.m[:yS]) # extract the load shedding
    x0 = value.(s.ex.m[:x0]) # extract the values of rhs parameters in the subproblem (eg, capacity)
    push!(s.q,Q__type(ϕ,c0,yG,yI,yO,yL,yS,x0)) # add new extended exact solution

    return s 
end

function update_m_!(s::S2_type)::S2_type

    # function to add a new exact solution to S2_type.m (exact solutions stored)
    # for given S2_type (Benders subproblem) -> algorithm 2 (Adapt_Bend)

    s.m.n += 1 # update number of solution stored (+1)
    s.m.θ[s.m.n] = s.temp.θl[s.temp.i] # store optimal solution θᵢ
    s.m.λ[s.m.n,:] .= s.temp.λl[:,s.temp.i] # store subgradient λᵢ wrt xᵢ
    s.m.ϕ[s.m.n,:] .= s.temp.ϕu[:,s.temp.i] # store subgradient ϕᵢ wrt cᵢ
    s.m.x[s.m.n,:] .= s.temp.x[:,s.temp.i] # store xᵢ
    s.m.c[s.m.n,:] .= s.temp.c[:,s.temp.i] # store cᵢ

    return s
end

function update_lb!(s::S2_type)::S2_type

    # function to add a new exact solution to S2_type.olb (oracle for lower bound)
    # for given S2_type (Benders subproblem) -> algorithm 2 (Adapt_Bend)

    s.olb.cons[s.m.n] = @constraint(s.olb.m, s.olb.m[:ϕ] + s.olb.m[:γ]'*s.m.c[s.m.n,:] >= 0. ) # add constraint (ϕ + γᵀcₙ >= 0) to the oracle
    s.olb.help[1:s.m.n] .= s.m.θ[1:s.m.n] .- sum!(ones(s.m.n),s.m.λ[1:s.m.n,:].*s.m.x[1:s.m.n,:]) # compute vector h used when solving the lb oracle (hₛ = θₛ - λₛᵀxₛ, ∀s=1,..,n)
    s.olb.n += 1 # update number of solution stored (+1)

    return s
end

function update_ub!(s::S2_type)::S2_type

    # function to add a new exact solution in S2_type.olb (oracle for upper bound)
    # for given S2_type (Benders subproblem) -> algorithm 2 (Adapt_Bend)

    s.oub.cons[s.m.n] = @constraint(s.oub.m, s.oub.m[:ϕ] + s.oub.m[:γ]'*s.m.x[s.m.n,:] <= 0. ) # add constraint (ϕ + γᵀxₙ <= 0) to the oracle
    s.oub.help[1:s.m.n] .= s.m.θ[1:s.m.n] .- sum!(ones(s.m.n),s.m.ϕ[1:s.m.n,:].*s.m.c[1:s.m.n,:]) # compute vector h used when solving the ub oracle (hₛ = θₛ - ϕₛᵀcₛ, ∀s)
    s.oub.n += 1 # update number of solution stored (+1)

    return s
end

function update_s!(s::S2_type)::S2_type

    # function to add a new exact solution to S2_type.m (exact solutions stored), S2_type.olb (oracle for lower bound), and S2_type.olb (oracle for upper bound)
    # for given S2_type (Benders subproblem) -> algorithm 2 (Adapt_Bend)

    update_m_!(s) # store new exact solution in S2_type.m
    update_lb!(s) # store new exact solution in lower bound oracle
    update_ub!(s) # store new exact solution in upper bound oracle

    return s
end

function run_lb!(s::S2_type)::S2_type

    # function to solve S2_type.olb (oracle for lower bound) and store the solutions in S2_type.temp (temporary data)
    # for given S2_type (Benders subproblem) -> algorithm 2 (Adapt_Bend)

    s.olb.objf = s.olb.vars'*vcat([1.],s.temp.c[:,s.temp.i]) # compute objective function for subproblem i (ϕ + γᵀcᵢ)
    @objective(s.olb.m, Min, s.olb.objf) # set objective function
    set_normalized_rhs.(s.olb.cons[1:s.olb.n], s.olb.help[1:s.olb.n].+ s.m.λ[1:s.olb.n,:]*s.temp.x[:,s.temp.i]) # update the rhs of the constraints (ϕ + γᵀcₛ >= hₛ + λₛᵀxᵢ, ∀s)
    optimize!(s.olb.m) # solve the JuMP model to optimality
    s.temp.θl[s.temp.i] = objective_value(s.olb.m) # store the valid lower bound θᵢ
    s.temp.λl[:,s.temp.i] .= s.m.λ[1:s.m.n,:]'*dual.(s.olb.cons[1:s.olb.n]) # store the valid subgradient λᵢ wrt xᵢ

    return s
end

function run_ub!(s::S2_type)::S2_type

    # function to solve S2_type.oub (oracle for upper bound) and store the solutions in S2_type.temp (temporary data)
    # for given S2_type (Benders subproblem) -> algorithm 2 (Adapt_Bend)

    s.oub.objf = s.oub.vars'*vcat([1.],s.temp.x[:,s.temp.i]) # compute objective function for subproblem i (ϕ + γᵀxᵢ)
    @objective(s.oub.m, Max, s.oub.objf) # set objective function
    set_normalized_rhs.(s.oub.cons[1:s.oub.n], s.oub.help[1:s.oub.n].+ s.m.ϕ[1:s.oub.n,:]*s.temp.c[:,s.temp.i])  # update the rhs of the constraints (ϕ + γᵀxₛ <= hₛ + ϕₛᵀcᵢ, ∀s)
    optimize!(s.oub.m) # solve the JuMP model to optimality
    s.temp.θu[s.temp.i] = objective_value(s.oub.m) # store the valid upper bound θᵢ
    s.temp.ϕu[:,s.temp.i] .= s.m.ϕ[1:s.m.n,:]'*dual.(s.oub.cons[1:s.oub.n]) # store the valid subgradient ϕᵢ wrt cᵢ

    return s
end

function run_oracles!(s::S2_type)::S2_type

    # function to call S2_type.olb (oracle for lower bound) and S2_type.oub (oracle for upper bound) to generate valid bounds on subproblem i
    # for given S2_type (Benders subproblem) -> algorithm 2 (Adapt_Bend)

    run_lb!(s) # call the oracle for lower bound at (xᵢ,cᵢ)
    run_ub!(s) # call the oracle for upper bound at (xᵢ,cᵢ)

    return s
end

# */ --------------------------------------------------------------------------------------------- /* #

function step_a!(b::B1_type,s::S1_type)::Tuple{B1_type,S1_type}

    # function to perform 'step_a' of the Benders algorithm (solve the relaxed master problem and store decisions {xᵢ, ∀i})
    # for given B1_type (Benders master problem) and S1_type (Benders subproblem) -> algorithm 1 (Stand_Bend) 

    optimize!(b.rmp.m) # solve the relaxed master problem to optimality
    b.temp.β .= value.(b.rmp.m[:β]) # store optimal values of {βᵢ, ∀i}
    b.temp.x .= value.(b.rmp.m[:x]) # store opitmal decisions {xᵢ, ∀i}
    s.temp.x .= b.temp.x # send decisions {xᵢ, ∀i} to the Benders subproblem 

    return b,s
end

function step_a!(b::B3_type,s::S3_type)::Tuple{B3_type,S3_type}

    # function to perform 'step_a' of the Benders algorithm (solve the relaxed master problem and store decisions {xᵢ, ∀i})
    # for given B3_type (Benders master problem) and S3_type (Benders subproblem) -> algorithm 3 (Zaker_Bend)

    optimize!(b.rmp.m) # solve the relaxed master problem to optimality
    b.temp.β .= value.(b.rmp.m[:β]) # store optimal values of {βᵢ, ∀i}
    b.temp.x .= value.(b.rmp.m[:x]) # store opitmal decisions {xᵢ, ∀i}
    s.temp.x .= b.temp.x # send decisions {xᵢ, ∀i} to the Benders subproblem 

    return b,s
end

function step_a!(b::B2_type,s::S2_type)::Tuple{B2_type,S2_type}

    # function to perform 'step_a' of the Benders algorithm (solve the relaxed master problem and store decisions {xᵢ, ∀i})
    # for given B2_type (Benders master problem) and S2_type (Benders subproblem) -> algorithm 2 (Adapt_Bend)

    optimize!(b.rmp.m) # solve the relaxed master problem to optimality
    b.temp.β .= value.(b.rmp.m[:β]) # store optimal values of {βᵢ, ∀i}
    b.temp.x .= value.(b.rmp.m[:x]) # store opitmal decisions {xᵢ, ∀i}
    s.temp.x .= b.temp.x # send decisions {xᵢ, ∀i} to the Benders subproblem 

    return b,s
end

function step_b!(b::B1_type)::B1_type

    # function to perform 'step_b' of the Benders algorithm (compute lower bound on optimal objective)
    # for given B1_type (Benders master problem) -> algorithm 1 (Stand_Bend) 

    b.hist.L[b.hist.k] = exp10(6)*objective_value(b.rmp.m) # compute and store the lower bound Lₖ on the optimal objective @ iter k

    return b
end

function step_b!(b::B3_type)::B3_type

    # function to perform 'step_b' of the Benders algorithm (compute lower bound on optimal objective)
    # for given B3_type (Benders master problem) -> algorithm 3 (Zaker_Bend)
    
    b.hist.L[b.hist.k] = exp10(6)*objective_value(b.rmp.m) # compute and store the lower bound Lₖ on the optimal objective @ iter k
    
    return b
end

function step_b!(b::B2_type)::B2_type

    # function to perform 'step_b' of the Benders algorithm (compute lower bound on optimal objective)
    # for given B2_type (Benders master problem) -> algorithm 2 (Adapt_Bend)

    b.hist.L[b.hist.k] = exp10(6)*objective_value(b.rmp.m) # compute and store the lower bound Lₖ on the optimal objective @ iter k

    return b
end

function step_c!(b::B1_type,s::S1_type)::Tuple{B1_type,S1_type}

    # function to perform 'step_c' of the Benders algorithm (solve each subproblem exactly and store optimal solutions)
    # for given B1_type (Benders master problem) and S1_type (Benders subproblem) -> algorithm 1 (Stand_Bend) 

    for s.temp.i in 1:b.data.unc.ni
        solv_exact!(s) # solve subproblem i exactly and store solutions
    end
    b.temp.θ .= s.temp.θ # send optimal objectives {θᵢ, ∀i} to the Benders master problem      
    b.temp.λ .= s.temp.λ # send subgradients {λᵢ, ∀i} to the Benders master problem     
    
    return b,s
end

function step_c!(b::B3_type,s::S3_type)::Tuple{B3_type,S3_type}

    # function to perform 'step_c' of the Benders algorithm (solve each subproblem to tolerance δ and store optimal solutions)
    # for given B3_type (Benders master problem) and S3_type (Benders subproblem) -> algorithm 3 (Zaker_Bend)

    δ = max(exp10(-8),b.data.q^(-b.hist.k+1.)) # compute optimality tolerance δ = q^(-k+1), where q is the damping factor and k is the current iter number
    JuMP.set_optimizer_attribute(s.ex.m,"BarConvTol",δ) # set optimality tolerance δ
    for s.temp.i in 1:b.data.unc.ni
        solv_exact!(s) # solve subproblem i up to tolerance δ and store solutions
    end
    b.temp.θp .= s.temp.θp # send primal objectives {θᵢ, ∀i} to the Benders master problem       
    b.temp.θd .= s.temp.θd # send dual objectives {θᵢ, ∀i} to the Benders master problem     
    b.temp.λ  .= s.temp.λ # send subgradients {λᵢ, ∀i} to the Benders master problem
    
    return b,s
end

function step_c!(b::B2_type,s::S2_type)::Tuple{B2_type,S2_type}

    # function to perform 'step_c' of the Benders algorithm (solve w subproblem exactly and store optimal solutions)
    # for given B2_type (Benders master problem) and S2_type (Benders subproblem) -> algorithm 2 (Adapt_Bend)

    set_Iex!(b) # select the w subproblems to solve exactly
    for s.temp.i in b.temp.Ie
        solv_exact!(s) # solve subproblem i exactly and store solutions
        update_s!(s) # send the new exact solution to the adaptive oracles
        extended_sol!(s) # store the new extended exact solution
    end

    return b,s
end

function step_d!()

end

function step_d!(b::B2_type,s::S2_type)::Tuple{B2_type,S2_type}

    # function to perform 'step_d' of the Benders algorithm (solve the adaptive oracles to generate inexact but valid solutions)
    # for given B2_type (Benders master problem) and S2_type (Benders subproblem) -> algorithm 2 (Adapt_Bend)

    for s.temp.i in 1:b.data.unc.ni
        run_oracles!(s) # solve the adaptive oracles for subproblem i and store the solutions
    end
    b.temp.θl .= s.temp.θl # send valid lower bound {θᵢ, ∀i} to the Benders master problem 
    b.temp.θu .= s.temp.θu # send valid upper bound {θᵢ, ∀i} to the Benders master problem
    b.temp.λl .= s.temp.λl # send valid subgradients {λᵢ, ∀i} to the Benders master problem
    b.temp.ϕu .= s.temp.ϕu # send valid subgradients {ϕᵢ, ∀i} to the Benders master problem

    return b,s
end

function step_e!(b::B1_type)::B1_type

    # function to perform 'step_e' of the Benders algorithm (compute upper bound on optimal objective)
    # for given B1_type (Benders master problem) -> algorithm 1 (Stand_Bend) 

    if (b.hist.k==1)
        b.hist.U[b.hist.k] = exp10(6)*(value(b.rmp.m[:f])+b.data.mp.κ*b.data.mp.π'*b.temp.θ) # compute and store the upper bound Uₖ on the optimal objective @ iter k = 1
    else
        b.hist.U[b.hist.k] = min(b.hist.U[b.hist.k-1],exp10(6)*(value(b.rmp.m[:f])+b.data.mp.κ*b.data.mp.π'*b.temp.θ)) # compute and store the upper bound Uₖ on the optimal objective @ iter k > 1
    end

    return b
end

function step_e!(b::B3_type)::B3_type

    # function to perform 'step_e' of the Benders algorithm (compute upper bound on optimal objective)
    # for given B3_type (Benders master problem) -> algorithm 3 (Zaker_Bend)

    if (b.hist.k==1)
        b.hist.U[b.hist.k] = exp10(6)*(value(b.rmp.m[:f])+b.data.mp.κ*b.data.mp.π'*b.temp.θp) # compute and store the upper bound Uₖ on the optimal objective @ iter k = 1
    else
        b.hist.U[b.hist.k] = min(b.hist.U[b.hist.k-1],exp10(6)*(value(b.rmp.m[:f])+b.data.mp.κ*b.data.mp.π'*b.temp.θp)) # compute and store the upper bound Uₖ on the optimal objective @ iter k > 1
    end

    return b
end

function step_e!(b::B2_type)::B2_type

    # function to perform 'step_e' of the Benders algorithm (compute upper bound on optimal objective)
    # for given B2_type (Benders master problem) -> algorithm 2 (Adapt_Bend)

    if (b.hist.k==1)
        b.hist.U[b.hist.k] = exp10(6)*(value(b.rmp.m[:f])+b.data.mp.κ*b.data.mp.π'*b.temp.θu) # compute and store the upper bound Uₖ on the optimal objective @ iter k = 1
    else
        b.hist.U[b.hist.k] = min(b.hist.U[b.hist.k-1], exp10(6)*(value(b.rmp.m[:f])+(b.data.mp.κ*b.data.mp.π'*b.temp.θu))) # compute and store the upper bound Uₖ on the optimal objective @ iter k > 1
    end

    return b
end

function step_f!(b::B1_type)::B1_type

    # function to perform 'step_f' of the Benders algorithm (compute optimality gap Δ and update the relaxed master problem)
    # for given B1_type (Benders master problem) -> algorithm 1 (Stand_Bend) 

    b.temp.Δ = (b.hist.U[b.hist.k]-b.hist.L[b.hist.k])/b.hist.U[b.hist.k]*100 # compute optimality gap Δ = (Uₖ-Lₖ)/Uₖ @ iter k
    for i in 1:b.data.unc.ni
        @constraint(b.rmp.m,b.rmp.m[:β][i] >= b.temp.θ[i]+b.temp.λ[:,i]'*(b.rmp.m[:x][:,i].-b.temp.x[:,i])) # add constraints βᵢ >= θᵢ + λᵢᵀ(𝑥ᵢ-xᵢ) to the relaxed master problem
    end
    
    return b
end

function step_f!(b::B3_type)::B3_type

    # function to perform 'step_f' of the Benders algorithm (compute optimality gap Δ and update the relaxed master problem)
    # for given B3_type (Benders master problem) -> algorithm 3 (Zaker_Bend)

    b.temp.Δ = (b.hist.U[b.hist.k]-b.hist.L[b.hist.k])/b.hist.U[b.hist.k]*100 # compute optimality gap Δ = (Uₖ-Lₖ)/Uₖ @ iter k
    for i in 1:b.data.unc.ni
        @constraint(b.rmp.m,b.rmp.m[:β][i] >= b.temp.θd[i]+b.temp.λ[:,i]'*(b.rmp.m[:x][:,i].-b.temp.x[:,i])) # add constraints βᵢ >= θᵢ + λᵢᵀ(𝑥ᵢ-xᵢ) to the relaxed master problem
    end
    
    return b
end

function step_f!(b::B2_type)::B2_type

    # function to perform 'step_f' of the Benders algorithm (compute optimality gap Δ and update the relaxed master problem)
    # for given B2_type (Benders master problem) -> algorithm 2 (Adapt_Bend)

    b.temp.Δ = (b.hist.U[b.hist.k]-b.hist.L[b.hist.k])/b.hist.U[b.hist.k]*100 # compute optimality gap Δ = (Uₖ-Lₖ)/Uₖ @ iter k
    for i in 1:b.data.unc.ni
        @constraint(b.rmp.m,b.rmp.m[:β][i] >= b.temp.θl[i]+b.temp.λl[:,i]'*(b.rmp.m[:x][:,i].-b.temp.x[:,i])) # add constraints βᵢ >= θᵢ + λᵢᵀ(𝑥ᵢ-xᵢ) to the relaxed master problem
    end
    
    return b
end

# */ --------------------------------------------------------------------------------------------- /* #

function step_0!(b::B2_type,s::S2_type)::Tuple{B2_type,S2_type}

    # function to perform 'step_0' of the Benders algorithm (solve subproblem exactly at special point (x_min,c_min) to make the adaptive oracles always feasible)
    # for given B2_type (Benders master problem) and S2_type (Benders subproblem) -> algorithm 2 (Adapt_Bend)

    s.temp.i = 1 # set i = 1
    s.temp.x[:,1] .= round.(vcat(minimum(b.data.mp.xh,dims=2)[:]*.99,minimum(b.data.unc.h,dims=1)[:].*[.99,1.01]); digits=4) # set x₁ = x_min
    s.temp.c[:,1] .= round.(minimum(b.data.unc.c,dims=1)[:]*.99; digits=4) # set c₁ = c_min
    solv_exact!(s) # solve subproblem 1 exactly (at point (x_min,c_min)) and store solutions
    update_s!(s) # send the new exact solution to the adaptive oracles
    s.temp.x .= vcat(b.data.mp.xh,b.data.unc.h') # set decisions {xᵢ, ∀i} on historial data
    s.temp.c .= b.data.unc.c' # set {cᵢ, ∀i}
    for s.temp.i in 1:b.data.unc.ni
        run_oracles!(s) # solve the adaptive oracles for subproblem i and store the solutions 
    end
    b.temp.θl .= s.temp.θl # send valid lower bound {θᵢ, ∀i} to the Benders master problem 
    b.temp.θu .= s.temp.θu # send valid upper bound {θᵢ, ∀i} to the Benders master problem
    b.temp.λl .= s.temp.λl # send valid subgradients {λᵢ, ∀i} to the Benders master problem
    b.temp.ϕu .= s.temp.ϕu # send valid subgradients {ϕᵢ, ∀i} to the Benders master problem

    return b,s
end

function step!(b::B1_type,s::S1_type)::Tuple{B1_type,S1_type}

    # function to perform a Benders 'step'
    # for given B1_type (Benders master problem) and S1_type (Benders subproblem) -> algorithm 1 (Stand_Bend) 

    b.hist.k += 1                                # update iter k number
    b.hist.T[b.hist.k,1] = @elapsed step_a!(b,s) # solve the relaxed master problem and store decisions {xᵢ, ∀i}
    b.hist.T[b.hist.k,2] = @elapsed step_b!(b)   # compute lower bound on optimal objective
    b.hist.T[b.hist.k,3] = @elapsed step_c!(b,s) # solve each subproblem exactly and store optimal solutions
    b.hist.T[b.hist.k,4] = @elapsed step_d!( )   #
    b.hist.T[b.hist.k,5] = @elapsed step_e!(b)   # compute upper bound on optimal objective
    b.hist.T[b.hist.k,6] = @elapsed step_f!(b)   # compute optimality gap Δ and update the relaxed master problem

    return b,s
end

function step!(b::B3_type,s::S3_type)::Tuple{B3_type,S3_type}

    # function to perform a Benders 'step'
    # for given B3_type (Benders master problem) and S3_type (Benders subproblem) -> algorithm 3 (Zaker_Bend)

    b.hist.k += 1                                # update iter k number
    b.hist.T[b.hist.k,1] = @elapsed step_a!(b,s) # solve the relaxed master problem and store decisions {xᵢ, ∀i}
    b.hist.T[b.hist.k,2] = @elapsed step_b!(b)   # compute lower bound on optimal objective
    b.hist.T[b.hist.k,3] = @elapsed step_c!(b,s) # solve each subproblem to tolerance δ and store optimal solutions
    b.hist.T[b.hist.k,4] = @elapsed step_d!( )   #
    b.hist.T[b.hist.k,5] = @elapsed step_e!(b)   # compute upper bound on optimal objective
    b.hist.T[b.hist.k,6] = @elapsed step_f!(b)   # compute optimality gap Δ and update the relaxed master problem

    return b,s
end

function step!(b::B2_type,s::S2_type)::Tuple{B2_type,S2_type}

    # function to perform a Benders 'step'
    # for given B2_type (Benders master problem) and S2_type (Benders subproblem) -> algorithm 2 (Adapt_Bend)

    b.hist.k += 1                                # update iter k number
    b.hist.T[b.hist.k,1] = @elapsed step_a!(b,s) # solve the relaxed master problem and store decisions {xᵢ, ∀i}
    b.hist.T[b.hist.k,2] = @elapsed step_b!(b)   # compute lower bound on optimal objective
    b.hist.T[b.hist.k,3] = @elapsed step_c!(b,s) # solve w subproblem exactly and store optimal solutions
    b.hist.T[b.hist.k,4] = @elapsed step_d!(b,s) # solve the adaptive oracles to generate inexact but valid solutions
    b.hist.T[b.hist.k,5] = @elapsed step_e!(b)   # compute upper bound on optimal objective
    b.hist.T[b.hist.k,6] = @elapsed step_f!(b)   # compute optimality gap Δ and update the relaxed master problem

    return b,s
end

# */ --------------------------------------------------------------------------------------------- /* #

function solve_deterministic!(b::JuMP.Model,c::Int64;print_output::Int64=1,save_results::Int64=1)::JuMP.Model

    # function to solve the deterministic equivalent of the stochastic investment planning problem
    # for given Model (master problem) and c (case)

    print_output == 1 ? print0(b,c) : nothing
    println(" Sending JuMP model to the Optimizer")
    println(" ")
    print_output == 1 ? nothing : JuMP.set_optimizer_attribute(b,"OutputFlag",0)
    optimize!(b)
    print_output == 1 ? print_summary_a(b,c) : nothing
    
    return b
end

function solve_Benders!(b::B1_type,s::S1_type;print_output::Int64=1)::Tuple{B1_type,S1_type}

    # function to run Stand_Bend algorithm to solve the investment planning problem
    # for given B1_type (Benders master problem) and S1_type (Benders subproblem) -> algorithm 1 (Stand_Bend) 
    
    print_output == 1 ? print0(b) : nothing
    while (b.hist.k < b.data.J) # check the iteration number k is lower than the limit J 
        step!(b,s) # perform a Benders step
        print_output == 1 ? try print_info(b) catch; end : nothing
        (b.temp.Δ <= b.data.δ) ? break : nothing # check if the optimality gap Δₖ is lower than the targer tolerance δ 
    end
    print_output == 1 ? print_summary(b) : nothing
    
    return b,s
end

function solve_Benders!(b::B3_type,s::S3_type;print_output::Int64=1)::Tuple{B3_type,S3_type}

    # function to run Zaker_Bend algorithm to solve the investment planning problem
    # for given B3_type (Benders master problem) and S3_type (Benders subproblem) -> algorithm 3 (Zaker_Bend) 
    
    print_output == 1 ? print0(b) : nothing
    while (b.hist.k < b.data.J) # check the iteration number k is lower than the limit J 
        step!(b,s) # perform a Benders step
        print_output == 1 ? try print_info(b) catch; end : nothing
        (b.temp.Δ <= b.data.δ) ? break : nothing # check if the optimality gap Δₖ is lower than the targer tolerance δ
    end
    print_output == 1 ? print_summary(b) : nothing
    
    return b,s
end

function solve_Benders!(b::B2_type,s::S2_type;print_output::Int64=1)::Tuple{B2_type,S2_type}
    
    print_output == 1 ? print0(b) : nothing
    step_0!(b,s) # perform step0 to solve the special point (x_min,c_min)
    while (b.hist.k < b.data.J) # check the iteration number k is lower than the limit J
        step!(b,s) # perform a Benders step
        print_output == 1 ? try print_info(b) catch; end : nothing
        (b.temp.Δ <= b.data.δ) ? break : nothing # check if the optimality gap Δₖ is lower than the targer tolerance δ
    end
    print_output == 1 ? print_summary(b) : nothing
    
    return b,s
end

# */ --------------------------------------------------------------------------------------------- /* #