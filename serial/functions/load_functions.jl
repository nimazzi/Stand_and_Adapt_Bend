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
    @variable(m, f, container=Array)
    @variable(m, x0[ms.P,ms.I0] >= .0, container=Array)
    @variable(m, x[xr,ms.I], container=Array)
    @variable(m, β[ms.I] >= .0, container=Array)
    # */ -- subproblems variables ---------------------------------- /* #
    @variable(m, ϕ[ms.I,cr], container=Array)
    @variable(m, c0[ms.I], container=Array) 
    @variable(m, yG[ms.I,ps.G,ps.S,ps.H] >= 0, container=Array)
    @variable(m, yI[ms.I,ps.B,ps.S,ps.H] >= 0, container=Array) 
    @variable(m, yO[ms.I,ps.B,ps.S,ps.H] >= 0, container=Array) 
    @variable(m, yL[ms.I,ps.B,ps.S,ps.H] >= 0, container=Array) 
    @variable(m, yS[ms.I,ps.S,ps.H] >= 0, container=Array) 
    @variable(m, x1[ms.I,xr], container=Array) 
    # */ ----------------------------------------------------------- /* #
    
    # */ -- obj function ------------------------------------------- /* #
    @objective(m, Min, f + mp.κ*sum(mp.π[i]*β[i] for i in ms.I) )
    # */ ----------------------------------------------------------- /* #
    
    # */ -- master constraints ------------------------------------- /* #
    @constraint(m, cm01, f >= exp10(-6)*sum(mp.π0[i0]*sum(mp.ci[p,i0]*x0[p,i0] for p in ms.P) for i0 in ms.I0) + exp10(-6)*mp.κ*sum(mp.π[i]*sum(mp.cf[p]*x[p,i] for p in ms.P) for i in ms.I) )
    @constraint(m, cm02[p=ms.P,i=ms.I], x[p,i] == mp.xh[p,i] + sum(x0[p,i0] for i0 in mp.map[i]) )
    @constraint(m, cm03[p=ms.P,i=ms.I], x[p,i] >= .0 )
    @constraint(m, cm04[p=ms.P,i=ms.I], x[p,i] <= mp.xm[p] )
    for i in ms.I for j in 1:unc.nh
        fix(x[unc.nx0+j,i],unc.h[i,j];force=true)
    end end
    # */ -- subproblems constraints -------------------------------- /* #
    @constraint(m, cs00[i=ms.I], β[i] >= c0[i] + unc.c[i,1]*ϕ[i,1] + unc.c[i,2]*ϕ[i,2])
    @constraint(m, cs01[i=ms.I], c0[i] == exp10(-6)*sum(sum(sum((pp.cvg[g]+pp.cfg[g]/pp.ηg[g])*yG[i,g,s,h] for g in ps.G0) + pp.cvg[ps.gn]*yG[i,ps.gn,s,h] + pp.cs*yS[i,s,h] for h in ps.H)*pp.α for s in ps.S))
    @constraint(m, cs02[i=ms.I], ϕ[i,1] == exp10(-6)*sum(sum(sum(pp.eg[g]/pp.ηg[g]*yG[i,g,s,h] for g in ps.G) for h in ps.H)*pp.α for s in ps.S))
    @constraint(m, cs03[i=ms.I], ϕ[i,2] == exp10(-6)*sum(sum(1/pp.ηg[ps.gn]*yG[i,ps.gn,s,h] for h in ps.H)*pp.α for s in ps.S))
    @constraint(m, cs04[i=ms.I,b=ps.B,s=ps.S,h=ps.H], yL[i,b,s,h] - yL[i,b,s,circ(h,ps.H)] == pp.ηb[b]*yI[i,b,s,h] - yO[i,b,s,h])
    @constraint(m, cs05[i=ms.I,g=ps.G,s=ps.S,h=ps.H], yG[i,g,s,h] - yG[i,g,s,circ(h,ps.H)] <= pp.rg[g]*x1[i,g]) 
    @constraint(m, cs06[i=ms.I,g=ps.G,s=ps.S,h=ps.H], yG[i,g,s,circ(h,ps.H)] - yG[i,g,s,h] <= pp.rg[g]*x1[i,g]) 
    @constraint(m, cs07[i=ms.I,g=ps.G,s=ps.S,h=ps.H], yG[i,g,s,h] <= x1[i,g]) 
    @constraint(m, cs08[i=ms.I,b=ps.B,s=ps.S,h=ps.H], yI[i,b,s,h] <= pp.pb[b]*x1[i,b+ps.b0])
    @constraint(m, cs09[i=ms.I,b=ps.B,s=ps.S,h=ps.H], yO[i,b,s,h] <= pp.pb[b]*x1[i,b+ps.b0])
    @constraint(m, cs10[i=ms.I,b=ps.B,s=ps.S,h=ps.H], yL[i,b,s,h] <= x1[i,b]) 
    @constraint(m, cs11[i=ms.I], sum(yG[i,g,s,h]*pp.eg[g]/pp.ηg[g] for g in ps.G for s in ps.S for h in ps.H) <= pp.lco*x1[i,unc.nx0+1]) 
    @constraint(m, cs12[i=ms.I,s=ps.S,h=ps.H], sum(yG[i,g,s,h] for g in ps.G) + sum(yO[i,b,s,h]-yI[i,b,s,h] for b in ps.B) + yS[i,s,h] >= -x1[i,unc.nx0+2]*pp.pd[s,h] - sum(pp.pr[r,s,h]*x1[i,r+ps.r0] for r in ps.R))
    @constraint(m, cs13[i=ms.I,n=1:unc.nx0], x1[i,n]*exp10(-3) == x[n,i])
    @constraint(m, cs14[i=ms.I,n=1:unc.nh ], x1[i,unc.nx0+n]   == x[unc.nx0+n,i])
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
    @variable(m, f, container=Array)
    @variable(m, x0[ms.P,ms.I0] >= .0, container=Array)
    @variable(m, x[xr,ms.I], container=Array)
    @variable(m, β[ms.I] >= .0, container=Array)
    # */ ----------------------------------------------------------- /* #
    
    # */ -- obj function ------------------------------------------- /* #
    @objective(m, Min, f + mp.κ*sum(mp.π[i]*β[i] for i in ms.I) )
    # */ ----------------------------------------------------------- /* #
    
    # */ -- constraints -------------------------------------------- /* #
    @constraint(m, c01, f >= exp10(-6)*sum(mp.π0[i0]*sum(mp.ci[p,i0]*x0[p,i0] for p in ms.P) for i0 in ms.I0) + exp10(-6)*mp.κ*sum(mp.π[i]*sum(mp.cf[p]*x[p,i] for p in ms.P) for i in ms.I) )
    @constraint(m, c02[p=ms.P,i=ms.I], x[p,i] == mp.xh[p,i] + sum(x0[p,i0] for i0 in mp.map[i]) )
    @constraint(m, c03[p=ms.P,i=ms.I], x[p,i] >= .0 )
    @constraint(m, c04[p=ms.P,i=ms.I], x[p,i] <= mp.xm[p] )
    for i in ms.I for j in 1:unc.nh
        fix(x[unc.nx0+j,i],unc.h[i,j];force=true)
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
    @variable(m, ϕ[1:unc.nc], container=Array)
    @variable(m, c0, container=Array) 
    @variable(m, yG[ps.G,ps.S,ps.H] >= 0, container=Array)
    @variable(m, yI[ps.B,ps.S,ps.H] >= 0, container=Array) 
    @variable(m, yO[ps.B,ps.S,ps.H] >= 0, container=Array) 
    @variable(m, yL[ps.B,ps.S,ps.H] >= 0, container=Array) 
    @variable(m, yS[ps.S,ps.H] >= 0, container=Array) 
    @variable(m, x0[1:unc.nx], container=Array) 
    @variable(m, x[1:unc.nx], container=Array) 
    # */ ----------------------------------------------------------- /* #
    
    # */ -- obj function ------------------------------------------- /* #
    @objective(m, Min, c0 )
    # */ ----------------------------------------------------------- /* #
    
    # */ -- constraints -------------------------------------------- /* #
    @constraint(m, c01, c0 == exp10(-6)*sum(sum(sum((pp.cvg[g]+pp.cfg[g]/pp.ηg[g])*yG[g,s,h] for g in ps.G0) + pp.cvg[ps.gn]*yG[ps.gn,s,h] + pp.cs*yS[s,h] for h in ps.H)*pp.α for s in ps.S))
    @constraint(m, c02, ϕ[1] == exp10(-6)*sum(sum(sum(pp.eg[g]/pp.ηg[g]*yG[g,s,h] for g in ps.G) for h in ps.H)*pp.α for s in ps.S))
    @constraint(m, c03, ϕ[2] == exp10(-6)*sum(sum(1/pp.ηg[ps.gn]*yG[ps.gn,s,h] for h in ps.H)*pp.α for s in ps.S))
    @constraint(m, c04[b=ps.B,s=ps.S,h=ps.H], yL[b,s,h] - yL[b,s,circ(h,ps.H)] == pp.ηb[b]*yI[b,s,h] - yO[b,s,h])
    @constraint(m, c05[g=ps.G,s=ps.S,h=ps.H], yG[g,s,h] - yG[g,s,circ(h,ps.H)] <= pp.rg[g]*x0[g]) 
    @constraint(m, c06[g=ps.G,s=ps.S,h=ps.H], yG[g,s,circ(h,ps.H)] - yG[g,s,h] <= pp.rg[g]*x0[g]) 
    @constraint(m, c07[g=ps.G,s=ps.S,h=ps.H], yG[g,s,h] <= x0[g]) 
    @constraint(m, c08[b=ps.B,s=ps.S,h=ps.H], yI[b,s,h] <= pp.pb[b]*x0[b+ps.b0])
    @constraint(m, c09[b=ps.B,s=ps.S,h=ps.H], yO[b,s,h] <= pp.pb[b]*x0[b+ps.b0])
    @constraint(m, c10[b=ps.B,s=ps.S,h=ps.H], yL[b,s,h] <= x0[b]) 
    @constraint(m, c11, sum(yG[g,s,h]*pp.eg[g]/pp.ηg[g] for g in ps.G for s in ps.S for h in ps.H) <= pp.lco*x0[unc.nx0+1]) 
    @constraint(m, c12[s=ps.S,h=ps.H], sum(yG[g,s,h] for g in ps.G) + sum(yO[b,s,h]-yI[b,s,h] for b in ps.B) + yS[s,h] >= -x0[unc.nx0+2]*pp.pd[s,h] - sum(pp.pr[r,s,h]*x0[r+ps.r0] for r in ps.R))
    @constraint(m, c13[n=1:unc.nx0], x0[n]*exp10(-3) == x[n])
    @constraint(m, c14[n=1:unc.nh ], x0[unc.nx0+n]   == x[unc.nx0+n])
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
# */ --- functions to generate structs : on main core -------------------------------------------- /* #
# */ --------------------------------------------------------------------------------------------- /* #

function gen_kiidx(u::u_type,w::Int64)::Tuple{Array{BitArray{1},1},Array{Array{Int64,1},1}}

    if w == 1
        kidx = Array{BitArray{1},1}(undef,0)
        push!(kidx,zeros(Int64,u.ni).==0)
        iidx = Array{Array{Int64,1},1}(undef,0)
        push!(iidx,[i for i in 1:u.ni])
    else
        x = zeros(u.nh+u.nc,u.ni)
        for j in 1:u.nh
            uh = u.h[:,j]
            hmin,hmax = minimum(uh),maximum(uh)
            (hmax > hmin) ? x[j,:] .= (uh.-hmin)./(hmax-hmin) :  x[j,:] .= zeros(u.ni)
        end
        for j in 1:u.nc
            uc = u.c[:,j]
            cmin,cmax = minimum(uc),maximum(uc)
            (cmax > cmin) ? x[j+u.nh,:] .= (uc.-cmin)./(cmax-cmin) : x[j+u.nh,:] .= zeros(u.ni)
        end
        kass = kmeans(x,w).assignments
        kidx = Array{BitArray{1},1}(undef,0)
        for j in 1:w
            push!(kidx,kass.==j)
        end
        idx = [i for i in 1:u.ni]
        iidx = Array{Array{Int64,1},1}(undef,0)
        for j in 1:w
            push!(iidx,idx[kidx[j]])
        end
    end

    return kidx,iidx
end

function gen_D__type(ms::ms_type,mp::mp_type,ps::ps_type,pp::pp_type,u::u_type,c::Int64,a::Int64,q::Float64,w::Int64,j::Int64,e::Float64)::D__type
    
    k,i = gen_kiidx(u,w)

    return D__type(ms,mp,ps,pp,u,c,a,q,w,k,i,j,e)
end

function gen_R__type(d::D__type)::R__type
    
    m = JuMP.Model(JuMP.optimizer_with_attributes(()->Gurobi.Optimizer(gurobi_env)))
    JuMP.set_optimizer_attribute(m,"OutputFlag",0)
    RMP!(m,d.ms,d.mp,d.unc)
    v = Array{JuMP.VariableRef,2}(undef,d.unc.nx+1,d.unc.ni)
    for i in 1:d.unc.ni v[:,i] .= vcat([m[:β][i]],m[:x][:,i]) end
    c = Array{JuMP.ConstraintRef,2}(undef,d.J,d.unc.ni)
    
    return R__type(m,v,c)
end

function gen_T1_type(d::D__type)::T1_type
    
    x = zeros(d.unc.nx,d.unc.ni)
    c = convert(Array{Float64,2},d.unc.c')
    β = zeros(d.unc.ni)
    θ = zeros(d.unc.ni)
    λ = zeros(d.unc.nx,d.unc.ni)
    
    return T1_type(x,c,β,θ,λ,0.)
end

function gen_T3_type(d::D__type)::T3_type
    
    x  = zeros(d.unc.nx,d.unc.ni)
    c  = convert(Array{Float64,2},d.unc.c')
    β  = zeros(d.unc.ni)
    θp = zeros(d.unc.ni)
    θd = zeros(d.unc.ni)
    λ  = zeros(d.unc.nx,d.unc.ni)
    
    return T3_type(x,c,β,θp,θd,λ,0.)
end

function gen_T2_type(d::D__type)::T2_type
    
    x  = zeros(d.unc.nx,d.unc.ni)
    c  = convert(Array{Float64,2},d.unc.c')
    β  = zeros(d.unc.ni)
    θl = zeros(d.unc.ni)
    λl = zeros(d.unc.nx,d.unc.ni)
    θu = zeros(d.unc.ni)
    ϕu = zeros(d.unc.nc,d.unc.ni)
    Ie = zeros(Int64,d.w)
    
    return T2_type(x,c,β,θl,λl,θu,ϕu,Ie,0.)
end

function gen_H__type(d::D__type)::H__type
    
    L = zeros(d.J)
    U = zeros(d.J)
    t = zeros(d.J)
    T = zeros(d.J,6)
    
    return H__type(L,U,t,T,0)
end

function gen_M__type(d::D__type)::M__type

    θ = zeros(nprocs())
    λ = zeros(nprocs(),d.unc.nx)
    ϕ = zeros(nprocs(),d.unc.nc)
    x = zeros(nprocs(),d.unc.nx)
    c = zeros(nprocs(),d.unc.nc)

    return M__type(0,θ,λ,ϕ,x,c)
end

function gen_B0_type(c::Int64)::JuMP.Model
    
    ms,mp,ps,pp,unc = load_data(c)
    m = JuMP.Model(JuMP.optimizer_with_attributes(()->Gurobi.Optimizer(gurobi_env)))
    println(" Attempting to build JuMP model of deterministic equivalent")
    println(" ")
    MP!(m,ms,mp,ps,pp,unc)
    println(" JuMP Model build successfully!")
    println(" ")
    
    return m
end

function gen_B1_type(c::Int64,a::Int64,j::Int64,e::Float64)::B1_type
    
    ms,mp,ps,pp,unc = load_data(c)
    d = gen_D__type(ms,mp,ps,pp,unc,c,a,.0,1,j,e)
    r = gen_R__type(d)
    t = gen_T1_type(d)
    h = gen_H__type(d)
    
    return B1_type(r,t,h,d) 
end

function gen_B3_type(c::Int64,a::Int64,q::Float64,j::Int64,e::Float64)::B3_type
    
    ms,mp,ps,pp,unc = load_data(c)
    d = gen_D__type(ms,mp,ps,pp,unc,c,a,q,1,j,e)
    r = gen_R__type(d)
    t = gen_T3_type(d)
    h = gen_H__type(d)
    
    return B3_type(r,t,h,d) 
end

function gen_B2_type(c::Int64,a::Int64,w::Int64,j::Int64,e::Float64)::B2_type
    
    ms,mp,ps,pp,unc = load_data(c)
    d = gen_D__type(ms,mp,ps,pp,unc,c,a,.0,w,j,e)
    r = gen_R__type(d)
    t = gen_T2_type(d)
    h = gen_H__type(d)
    m = gen_M__type(d)
    
    return B2_type(r,t,h,m,d)
end

# */ --------------------------------------------------------------------------------------------- /* #
# */ --- functions to generate structs : everywhere ---------------------------------------------- /* #
# */ --------------------------------------------------------------------------------------------- /* #

function gen_d__type(d::D__type)::d__type

    s = d.ps
    p = d.pp
    u = d.unc
    a = d.algm
    j = d.J
    
    return d__type(s,p,u,a,j)
end

function gen_Lb_type(d::d__type)::O__type

    m = JuMP.Model(JuMP.optimizer_with_attributes(()->Gurobi.Optimizer(gurobi_env)))
    JuMP.set_optimizer_attribute(m,"OutputFlag",0)
    LBoracle!(m,d.unc,d.J+1)
    v = vcat(m[:ϕ],m[:γ])
    c = Array{JuMP.ConstraintRef,1}(undef,d.J+1)
    o = v'*ones(length(v))
    e = v'*ones(length(v))
    h = zeros(d.J+1)

    return O__type(m,c,v,o,e,h,0)
end

function gen_Ub_type(d::d__type)::O__type

    m = JuMP.Model(JuMP.optimizer_with_attributes(()->Gurobi.Optimizer(gurobi_env)))
    JuMP.set_optimizer_attribute(m,"OutputFlag",0)
    UBoracle!(m,d.unc,d.J+1)
    v = vcat(m[:ϕ],m[:γ])
    c = Array{JuMP.ConstraintRef,1}(undef,d.J+1)
    o = v'*ones(length(v))
    e = v'*ones(length(v))
    h = zeros(d.J+1)

    return O__type(m,c,v,o,e,h,0)
end

function gen_E__type(d::d__type)::E__type

    m = JuMP.Model(JuMP.optimizer_with_attributes(()->Gurobi.Optimizer(gurobi_env)))
    JuMP.set_optimizer_attribute(m,"OutputFlag",0)
    JuMP.set_optimizer_attribute(m,"Method",2)
    d.algm == 3 ? JuMP.set_optimizer_attribute(m,"Crossover",0) : nothing
    SP!(m,d.ps,d.pp,d.unc)
    v = vcat(m[:c0],m[:ϕ][:])
    o = v'*ones(length(v))

    return E__type(m,v,o)
end

function gen_M__type(d::d__type)::M__type

    θ = zeros(d.J+1)
    λ = zeros(d.J+1,d.unc.nx)
    ϕ = zeros(d.J+1,d.unc.nc)
    x = zeros(d.J+1,d.unc.nx)
    c = zeros(d.J+1,d.unc.nc)

    return M__type(0,θ,λ,ϕ,x,c)
end

function gen_t1_type(d::d__type)::t1_type

    x = zeros(d.unc.nx,d.unc.ni)
    c = convert(Array{Float64,2},d.unc.c')
    θ = zeros(d.unc.ni)
    λ = zeros(d.unc.nx,d.unc.ni)

   return t1_type(x,c,θ,λ,0)
end

function gen_t3_type(d::d__type)::t3_type

    x  = zeros(d.unc.nx,d.unc.ni)
    c  = convert(Array{Float64,2},d.unc.c')
    θp = zeros(d.unc.ni)
    θd = zeros(d.unc.ni)
    λ  = zeros(d.unc.nx,d.unc.ni)

   return t3_type(x,c,θp,θd,λ,0)
end

function gen_t2_type(d::d__type)::t2_type

    x  = zeros(d.unc.nx,d.unc.ni)
    c  = convert(Array{Float64,2},d.unc.c')
    θl = zeros(d.unc.ni)
    θu = zeros(d.unc.ni)
    λl = zeros(d.unc.nx,d.unc.ni)
    ϕu = zeros(d.unc.nc,d.unc.ni)

   return t2_type(x,c,θl,θu,λl,ϕu,0)
end

function gen_S1_type(D::D__type)::S1_type

    d = gen_d__type(D)
    e = gen_E__type(d)
    t = gen_t1_type(d)

    return S1_type(e,d,t)
end

function gen_S3_type(D::D__type)::S3_type

    d = gen_d__type(D)
    e = gen_E__type(d)
    t = gen_t3_type(d)

    return S3_type(e,d,t)
end

function gen_S2_type(D::D__type)::S2_type

    d = gen_d__type(D)
    e = gen_E__type(d)
    l = gen_Lb_type(d)
    u = gen_Ub_type(d)
    m = gen_M__type(d)
    t = gen_t2_type(d)

    return S2_type(e,l,u,m,d,t)
end

# */ --------------------------------------------------------------------------------------------- /* #
# */ --- functions to generate structs : on main core and everywhere ----------------------------- /* #
# */ --------------------------------------------------------------------------------------------- /* #

function generate_structures(c::Int64,a::Int64,q::Float64,w::Int64,j::Int64,e::Float64)

    if a == 0
        B = gen_B0_type(c)
        S = nothing
    end
    if a == 1
        B = gen_B1_type(c,a,j,e)
        S = gen_S1_type(B.data)
    end
    if a == 2
        B = gen_B2_type(c,a,w,j,e)
        S = gen_S2_type(B.data)
    end
    if a == 3
        B = gen_B3_type(c,a,q,j,e)
        S = gen_S3_type(B.data)
    end
    
    return B,S
end

# */ --------------------------------------------------------------------------------------------- /* #
# */ --- lower level functions : on main core ---------------------------------------------------- /* #
# */ --------------------------------------------------------------------------------------------- /* #

function set_Iex!(b::B2_type)::B2_type

    for w in 1:b.data.w
        ie = findmax(b.data.mp.π[b.data.kidx[w]].*(b.temp.θu[b.data.kidx[w]].-b.temp.θl[b.data.kidx[w]]))[2]
        b.temp.Ie[w] = b.data.iidx[w][ie]
    end

    return b
end

# */ --------------------------------------------------------------------------------------------- /* #
# */ --- lower level functions : everywhere ------------------------------------------------------ /* #
# */ --------------------------------------------------------------------------------------------- /* #

function solv_exact!(s::S1_type)::S1_type

    s.ex.objf =  s.ex.vars'*vcat([1.],s.temp.c[:,s.temp.i])
    @objective(s.ex.m, Min, s.ex.objf)
    fix.(s.ex.m[:x],s.temp.x[:,s.temp.i];force=true)
    optimize!(s.ex.m)
    s.temp.θ[s.temp.i] = objective_value(s.ex.m)
    s.temp.λ[:,s.temp.i] .= dual.(FixRef.(s.ex.m[:x]))

    return s
end

function solv_exact!(s::S3_type)::S3_type

    s.ex.objf =  s.ex.vars'*vcat([1.],s.temp.c[:,s.temp.i])
    @objective(s.ex.m, Min, s.ex.objf)
    fix.(s.ex.m[:x],s.temp.x[:,s.temp.i];force=true)
    optimize!(s.ex.m)
    s.temp.θp[s.temp.i] = objective_value(s.ex.m)
    s.temp.λ[:,s.temp.i] .= dual.(FixRef.(s.ex.m[:x]))
    s.temp.θd[s.temp.i] = s.temp.λ[:,s.temp.i]'*value.(s.ex.m[:x])

    return s
end

function solv_exact!(s::S2_type)::S2_type

    s.ex.objf =  s.ex.vars'*vcat([1.],s.temp.c[:,s.temp.i])
    @objective(s.ex.m, Min, s.ex.objf)
    fix.(s.ex.m[:x],s.temp.x[:,s.temp.i];force=true)
    optimize!(s.ex.m)
    s.temp.θl[s.temp.i] = objective_value(s.ex.m)
    s.temp.λl[:,s.temp.i] .= dual.(FixRef.(s.ex.m[:x]))
    s.temp.ϕu[:,s.temp.i] .= value.(s.ex.m[:ϕ])

    return s
end

function update_m_!(s::S2_type)::S2_type

    s.m.n += 1
    s.m.θ[s.m.n] = s.temp.θl[s.temp.i]
    s.m.λ[s.m.n,:] .= s.temp.λl[:,s.temp.i]
    s.m.ϕ[s.m.n,:] .= s.temp.ϕu[:,s.temp.i]
    s.m.x[s.m.n,:] .= s.temp.x[:,s.temp.i]
    s.m.c[s.m.n,:] .= s.temp.c[:,s.temp.i]

    return s
end

function update_lb!(s::S2_type)::S2_type

    s.olb.cons[s.m.n] = @constraint(s.olb.m, s.olb.m[:ϕ] + s.olb.m[:γ]'*s.m.c[s.m.n,:] >= 0. )
    s.olb.help[1:s.m.n] .= s.m.θ[1:s.m.n] .- sum!(ones(s.m.n),s.m.λ[1:s.m.n,:].*s.m.x[1:s.m.n,:])
    s.olb.n += 1

    return s
end

function update_ub!(s::S2_type)::S2_type

    s.oub.cons[s.m.n] = @constraint(s.oub.m, s.oub.m[:ϕ] + s.oub.m[:γ]'*s.m.x[s.m.n,:] <= 0. )
    s.oub.help[1:s.m.n] .= s.m.θ[1:s.m.n] .- sum!(ones(s.m.n),s.m.ϕ[1:s.m.n,:].*s.m.c[1:s.m.n,:])
    s.oub.n += 1

    return s
end

function update_s!(s::S2_type)::S2_type

    update_m_!(s)
    update_lb!(s)
    update_ub!(s)

    return s
end

function run_lb!(s::S2_type)::S2_type

    s.olb.objf = s.olb.vars'*vcat([1.],s.temp.c[:,s.temp.i])
    @objective(s.olb.m, Min, s.olb.objf)
    set_normalized_rhs.(s.olb.cons[1:s.olb.n], s.olb.help[1:s.olb.n].+ s.m.λ[1:s.olb.n,:]*s.temp.x[:,s.temp.i])
    optimize!(s.olb.m)
    s.temp.θl[s.temp.i] = objective_value(s.olb.m)
    s.temp.λl[:,s.temp.i] .= s.m.λ[1:s.m.n,:]'*dual.(s.olb.cons[1:s.olb.n])

    return s
end

function run_ub!(s::S2_type)::S2_type

    s.oub.objf = s.oub.vars'*vcat([1.],s.temp.x[:,s.temp.i])
    @objective(s.oub.m, Max, s.oub.objf)
    set_normalized_rhs.(s.oub.cons[1:s.oub.n], s.oub.help[1:s.oub.n].+ s.m.ϕ[1:s.oub.n,:]*s.temp.c[:,s.temp.i])
    optimize!(s.oub.m)
    s.temp.θu[s.temp.i] = objective_value(s.oub.m)
    s.temp.ϕu[:,s.temp.i] .= s.m.ϕ[1:s.m.n,:]'*dual.(s.oub.cons[1:s.oub.n])

    return s
end

function run_oracles!(s::S2_type)::S2_type

    run_lb!(s)
    run_ub!(s)

    return s
end

# */ --------------------------------------------------------------------------------------------- /* #

function step_a!(b::B1_type,s::S1_type)::Tuple{B1_type,S1_type}

    optimize!(b.rmp.m)
    b.temp.β .= value.(b.rmp.m[:β])
    b.temp.x .= value.(b.rmp.m[:x])
    s.temp.x .= b.temp.x

    return b,s
end

function step_a!(b::B3_type,s::S3_type)::Tuple{B3_type,S3_type}

    optimize!(b.rmp.m)
    b.temp.β .= value.(b.rmp.m[:β])
    b.temp.x .= value.(b.rmp.m[:x])
    s.temp.x .= b.temp.x

    return b,s
end

function step_a!(b::B2_type,s::S2_type)::Tuple{B2_type,S2_type}

    optimize!(b.rmp.m)
    b.temp.β .= value.(b.rmp.m[:β])
    b.temp.x .= value.(b.rmp.m[:x])
    s.temp.x .= b.temp.x

    return b,s
end

function step_b!(b::B1_type)::B1_type

    b.hist.L[b.hist.k] = exp10(6)*objective_value(b.rmp.m)

    return b
end

function step_b!(b::B3_type)::B3_type
    
    b.hist.L[b.hist.k] = exp10(6)*objective_value(b.rmp.m)
    
    return b
end

function step_b!(b::B2_type)::B2_type

    b.hist.L[b.hist.k] = exp10(6)*objective_value(b.rmp.m)

    return b
end

function step_c!(b::B1_type,s::S1_type)::Tuple{B1_type,S1_type}

    for s.temp.i in 1:b.data.unc.ni
        solv_exact!(s)
    end
    b.temp.θ .= s.temp.θ      
    b.temp.λ .= s.temp.λ      
    
    return b,s
end

function step_c!(b::B3_type,s::S3_type)::Tuple{B3_type,S3_type}

    JuMP.set_optimizer_attribute(s.ex.m,"BarConvTol",max(exp10(-8),b.data.q^(-b.hist.k+1.)))
    for s.temp.i in 1:b.data.unc.ni
        solv_exact!(s)
    end
    b.temp.θp .= s.temp.θp      
    b.temp.θd .= s.temp.θd      
    b.temp.λ  .= s.temp.λ 
    
    return b,s
end

function step_c!(b::B2_type,s::S2_type)::Tuple{B2_type,S2_type}

    set_Iex!(b)
    for s.temp.i in b.temp.Ie
        solv_exact!(s)
        update_s!(s)
    end

    return b,s
end

function step_d!()

end

function step_d!(b::B2_type,s::S2_type)::Tuple{B2_type,S2_type}

    for s.temp.i in 1:b.data.unc.ni
        run_oracles!(s)
    end
    b.temp.θl .= s.temp.θl
    b.temp.θu .= s.temp.θu
    b.temp.λl .= s.temp.λl
    b.temp.ϕu .= s.temp.ϕu

    return b,s
end

function step_e!(b::B1_type)::B1_type

    if (b.hist.k==1)
        b.hist.U[b.hist.k] = exp10(6)*(value(b.rmp.m[:f])+b.data.mp.κ*b.data.mp.π'*b.temp.θ)
    else
        b.hist.U[b.hist.k] = min(b.hist.U[b.hist.k-1],exp10(6)*(value(b.rmp.m[:f])+b.data.mp.κ*b.data.mp.π'*b.temp.θ))
    end

    return b
end

function step_e!(b::B3_type)::B3_type

    if (b.hist.k==1)
        b.hist.U[b.hist.k] = exp10(6)*(value(b.rmp.m[:f])+b.data.mp.κ*b.data.mp.π'*b.temp.θp)
    else
        b.hist.U[b.hist.k] = min(b.hist.U[b.hist.k-1],exp10(6)*(value(b.rmp.m[:f])+b.data.mp.κ*b.data.mp.π'*b.temp.θp))
    end

    return b
end

function step_e!(b::B2_type)::B2_type

    if (b.hist.k==1)
        b.hist.U[b.hist.k] = exp10(6)*(value(b.rmp.m[:f])+b.data.mp.κ*b.data.mp.π'*b.temp.θu)
    else
        b.hist.U[b.hist.k] = min(b.hist.U[b.hist.k-1], exp10(6)*(value(b.rmp.m[:f])+(b.data.mp.κ*b.data.mp.π'*b.temp.θu)))
    end

    return b
end

function step_f!(b::B1_type)::B1_type

    b.temp.Δ = (b.hist.U[b.hist.k]-b.hist.L[b.hist.k])/b.hist.U[b.hist.k]*100
    for i in 1:b.data.unc.ni
        @constraint(b.rmp.m,b.rmp.m[:β][i] >= b.temp.θ[i]+b.temp.λ[:,i]'*(b.rmp.m[:x][:,i].-b.temp.x[:,i]))
    end
    
    return b
end

function step_f!(b::B3_type)::B3_type

    b.temp.Δ = (b.hist.U[b.hist.k]-b.hist.L[b.hist.k])/b.hist.U[b.hist.k]*100
    for i in 1:b.data.unc.ni
        @constraint(b.rmp.m,b.rmp.m[:β][i] >= b.temp.θd[i]+b.temp.λ[:,i]'*(b.rmp.m[:x][:,i].-b.temp.x[:,i]))
    end
    
    return b
end

function step_f!(b::B2_type)::B2_type

    b.temp.Δ = (b.hist.U[b.hist.k]-b.hist.L[b.hist.k])/b.hist.U[b.hist.k]*100
    for i in 1:b.data.unc.ni
        @constraint(b.rmp.m,b.rmp.m[:β][i] >= b.temp.θl[i]+b.temp.λl[:,i]'*(b.rmp.m[:x][:,i].-b.temp.x[:,i]))
    end
    
    return b
end

# */ --------------------------------------------------------------------------------------------- /* #

function step_0!(b::B2_type,s::S2_type)::Tuple{B2_type,S2_type}

    s.temp.i = 1
    s.temp.x[:,1] .= round.(vcat(minimum(b.data.mp.xh,dims=2)[:]*.99,minimum(b.data.unc.h,dims=1)[:].*[.99,1.01]); digits=4)
    s.temp.c[:,1] .= round.(minimum(b.data.unc.c,dims=1)[:]*.99; digits=4)
    solv_exact!(s)
    update_s!(s)
    s.temp.x .= vcat(b.data.mp.xh,b.data.unc.h')
    s.temp.c .= b.data.unc.c'
    for s.temp.i in 1:b.data.unc.ni
        run_oracles!(s)
    end
    b.temp.θl .= s.temp.θl
    b.temp.θu .= s.temp.θu
    b.temp.λl .= s.temp.λl
    b.temp.ϕu .= s.temp.ϕu

    return b,s
end

function step!(b::B1_type,s::S1_type)::Tuple{B1_type,S1_type}

    b.hist.k += 1
    b.hist.T[b.hist.k,1] = @elapsed step_a!(b,s)
    b.hist.T[b.hist.k,2] = @elapsed step_b!(b)
    b.hist.T[b.hist.k,3] = @elapsed step_c!(b,s)
    b.hist.T[b.hist.k,4] = @elapsed step_d!( )
    b.hist.T[b.hist.k,5] = @elapsed step_e!(b)
    b.hist.T[b.hist.k,6] = @elapsed step_f!(b)

    return b,s
end

function step!(b::B3_type,s::S3_type)::Tuple{B3_type,S3_type}

    b.hist.k += 1
    b.hist.T[b.hist.k,1] = @elapsed step_a!(b,s)
    b.hist.T[b.hist.k,2] = @elapsed step_b!(b)
    b.hist.T[b.hist.k,3] = @elapsed step_c!(b,s)
    b.hist.T[b.hist.k,4] = @elapsed step_d!( )
    b.hist.T[b.hist.k,5] = @elapsed step_e!(b)
    b.hist.T[b.hist.k,6] = @elapsed step_f!(b)

    return b,s
end

function step!(b::B2_type,s::S2_type)::Tuple{B2_type,S2_type}

    b.hist.k += 1
    b.hist.T[b.hist.k,1] = @elapsed step_a!(b,s)
    b.hist.T[b.hist.k,2] = @elapsed step_b!(b)
    b.hist.T[b.hist.k,3] = @elapsed step_c!(b,s)
    b.hist.T[b.hist.k,4] = @elapsed step_d!(b,s)
    b.hist.T[b.hist.k,5] = @elapsed step_e!(b)
    b.hist.T[b.hist.k,6] = @elapsed step_f!(b)

    return b,s
end

# */ --------------------------------------------------------------------------------------------- /* #

function solve_deterministic!(b::JuMP.Model,c::Int64;print_output::Int64=1,save_results::Int64=1)::JuMP.Model
    
    print_output == 1 ? print0(b,c) : nothing
    println(" Sending JuMP model to the Optimizer")
    println(" ")
    print_output == 1 ? nothing : JuMP.set_optimizer_attribute(b,"OutputFlag",0)
    optimize!(b)
    print_output == 1 ? print_summary_a(b,c) : nothing
    
    return b
end

function solve_Benders!(b::B1_type,s::S1_type;print_output::Int64=1)::Tuple{B1_type,S1_type}
    
    print_output == 1 ? print0(b) : nothing
    while (b.hist.k < b.data.J)
        step!(b,s)
        print_output == 1 ? try print_info(b) catch; end : nothing
        (b.temp.Δ <= b.data.δ) ? break : nothing
    end
    print_output == 1 ? print_summary(b) : nothing
    
    return b,s
end

function solve_Benders!(b::B3_type,s::S3_type;print_output::Int64=1)::Tuple{B3_type,S3_type}
    
    print_output == 1 ? print0(b) : nothing
    while (b.hist.k < b.data.J)
        step!(b,s)
        print_output == 1 ? try print_info(b) catch; end : nothing
        (b.temp.Δ <= b.data.δ) ? break : nothing
    end
    print_output == 1 ? print_summary(b) : nothing
    
    return b,s
end

function solve_Benders!(b::B2_type,s::S2_type;print_output::Int64=1)::Tuple{B2_type,S2_type}
    
    print_output == 1 ? print0(b) : nothing
    step_0!(b,s)
    while (b.hist.k < b.data.J)
        step!(b,s)
        print_output == 1 ? try print_info(b) catch; end : nothing
        (b.temp.Δ <= b.data.δ) ? break : nothing
    end
    print_output == 1 ? print_summary(b) : nothing
    
    return b,s
end

# */ --------------------------------------------------------------------------------------------- /* #