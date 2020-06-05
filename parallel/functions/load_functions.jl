# */ --------------------------------------------------------------------------------------------- /* #
# */ --- help functions -------------------------------------------------------------------------- /* #
# */ --------------------------------------------------------------------------------------------- /* #

@everywhere function circ(h::Int64,H::UnitRange{Int64})::Int64
    
    (h==minimum(H)) ? idx = maximum(H) : idx = h-1
    
    return idx
end

# */ --------------------------------------------------------------------------------------------- /* #
# */ --- optimization models --------------------------------------------------------------------- /* #
# */ --------------------------------------------------------------------------------------------- /* #

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

@everywhere function SP!(m::JuMP.Model,ps::ps_type,pp::pp_type,unc::u_type)::JuMP.Model
    
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

@everywhere function LBoracle!(m::JuMP.Model,unc::u_type)::JuMP.Model
    
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

@everywhere function UBoracle!(m::JuMP.Model,unc::u_type)::JuMP.Model
    
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

function gen_D__type(ms::ms_type,mp::mp_type,ps::ps_type,pp::pp_type,u::u_type,c::Int64,a::Int64,p::Int64,j::Int64,e::Float64)::D__type
    
    if nworkers()>1
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
        kass = kmeans(x,nworkers()).assignments
        kidx = Array{BitArray{1},1}(undef,0)
        for j in 1:nworkers()
            push!(kidx,kass.==j)
        end
        idx = [i for i in 1:u.ni]
        iidx = Array{Array{Int64,1},1}(undef,0)
        for j in 1:nworkers()
            push!(iidx,idx[kidx[j]])
        end
    else
        kidx = Array{BitArray{1},1}(undef,0)
        push!(kidx,zeros(Int64,u.ni).==0)
        iidx = Array{Array{Int64,1},1}(undef,0)
        push!(iidx,[i for i in 1:u.ni])
    end

    return D__type(ms,mp,ps,pp,u,c,a,p,kidx,iidx,j,e)
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
    re = Array{Tuple{Float64,Array{Float64,1}},1}(undef,d.unc.ni) 
    
    return T1_type(x,c,β,θ,λ,re,0.)
end

function gen_T2_type(d::D__type)::T2_type
    
    x  = zeros(d.unc.nx,d.unc.ni)
    c  = convert(Array{Float64,2},d.unc.c')
    β  = zeros(d.unc.ni)
    θl = zeros(d.unc.ni)
    λl = zeros(d.unc.nx,d.unc.ni)
    θu = zeros(d.unc.ni)
    ϕu = zeros(d.unc.nc,d.unc.ni)
    re = Array{Tuple{Float64,Array{Float64,1},Array{Float64,1}},1}(undef,nworkers()) 
    ro = Array{Tuple{Float64,Float64,Array{Float64,1},Array{Float64,1}},1}(undef,d.unc.ni) 
    Ie = zeros(Int64,0)
    
    return T2_type(x,c,β,θl,λl,θu,ϕu,re,ro,Ie,0.)
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

function gen_B1_type(c::Int64,a::Int64,p::Int64,j::Int64,e::Float64)::B1_type
    
    ms,mp,ps,pp,unc = load_data(c)
    d = gen_D__type(ms,mp,ps,pp,unc,c,a,p,j,e)
    r = gen_R__type(d)
    t = gen_T1_type(d)
    h = gen_H__type(d)
    
    return B1_type(r,t,h,d) 
end

function gen_B2_type(c::Int64,a::Int64,p::Int64,j::Int64,e::Float64)::B2_type
    
    ms,mp,ps,pp,unc = load_data(c)
    d = gen_D__type(ms,mp,ps,pp,unc,c,a,p,j,e)
    r = gen_R__type(d)
    t = gen_T2_type(d)
    h = gen_H__type(d)
    m = gen_M__type(d)
    
    return B2_type(r,t,h,m,d)
end

# */ --------------------------------------------------------------------------------------------- /* #
# */ --- functions to generate structs : everywhere ---------------------------------------------- /* #
# */ --------------------------------------------------------------------------------------------- /* #

@everywhere function gen_d__type(d::D__type)::d__type

    s = d.ps
    p = d.pp
    u = d.unc
    w = d.wrks
    j = d.J*nworkers()
    
    return d__type(s,p,u,w,j)
end

@everywhere function gen_Lb_type(d::d__type)::O__type

    m = JuMP.Model(JuMP.optimizer_with_attributes(()->Gurobi.Optimizer(gurobi_env)))
    JuMP.set_optimizer_attribute(m,"OutputFlag",0)
    JuMP.set_optimizer_attribute(m,"Threads",1)
    LBoracle!(m,d.unc)
    v = vcat(m[:ϕ],m[:γ])
    c = Array{JuMP.ConstraintRef,1}(undef,d.Jw+1)
    o = v'*ones(length(v))
    e = v'*ones(length(v))
    h = zeros(d.Jw+1)

    return O__type(m,c,v,o,e,h,0)
end

@everywhere function gen_Ub_type(d::d__type)::O__type

    m = JuMP.Model(JuMP.optimizer_with_attributes(()->Gurobi.Optimizer(gurobi_env)))
    JuMP.set_optimizer_attribute(m,"OutputFlag",0)
    JuMP.set_optimizer_attribute(m,"Threads",1)
    UBoracle!(m,d.unc)
    v = vcat(m[:ϕ],m[:γ])
    c = Array{JuMP.ConstraintRef,1}(undef,d.Jw+1)
    o = v'*ones(length(v))
    e = v'*ones(length(v))
    h = zeros(d.Jw+1)

    return O__type(m,c,v,o,e,h,0)
end

@everywhere function gen_E__type(d::d__type)::E__type

    m = JuMP.Model(JuMP.optimizer_with_attributes(()->Gurobi.Optimizer(gurobi_env)))
    JuMP.set_optimizer_attribute(m,"OutputFlag",0)
    JuMP.set_optimizer_attribute(m,"Threads",1)
    JuMP.set_optimizer_attribute(m,"Method",2)
    SP!(m,d.ps,d.pp,d.unc)
    v = vcat(m[:c0],m[:ϕ][:])
    o = v'*ones(length(v))

    return E__type(m,v,o)
end

@everywhere function gen_M__type(d::d__type)::M__type

    θ = zeros(d.Jw+1)
    λ = zeros(d.Jw+1,d.unc.nx)
    ϕ = zeros(d.Jw+1,d.unc.nc)
    x = zeros(d.Jw+1,d.unc.nx)
    c = zeros(d.Jw+1,d.unc.nc)

    return M__type(0,θ,λ,ϕ,x,c)
end

@everywhere function gen_t__type(d::d__type)::t__type

    x  = zeros(d.unc.nx,d.unc.ni)
    c  = convert(Array{Float64,2},d.unc.c')

   return t__type(x,c)
end

@everywhere function gen_S1_type(D::D__type)::S1_type

    d = gen_d__type(D)
    e = gen_E__type(d)
    t = gen_t__type(d)

    return S1_type(e,d,t)
end

@everywhere function gen_S2_type(D::D__type)::S2_type

    d = gen_d__type(D)
    e = gen_E__type(d)
    l = gen_Lb_type(d)
    u = gen_Ub_type(d)
    m = gen_M__type(d)
    t = gen_t__type(d)

    return S2_type(e,l,u,m,d,t)
end

# */ --------------------------------------------------------------------------------------------- /* #
# */ --- functions to generate structs : on main core and everywhere ----------------------------- /* #
# */ --------------------------------------------------------------------------------------------- /* #

function generate_structures(c::Int64,a::Int64,p::Int64,j::Int64,e::Float64)

    if a == 1
        B = gen_B1_type(c,a,p,j,e)
        sendto(procs(),d=B.data)
        @everywhere S = gen_S1_type(d)
    end
    if a == 2
        B = gen_B2_type(c,a,p,j,e)
        sendto(procs(),d=B.data)
        @everywhere S = gen_S2_type(d)
    end
    
    return B,S
end

# */ --------------------------------------------------------------------------------------------- /* #
# */ --- lower level functions : on main core ---------------------------------------------------- /* #
# */ --------------------------------------------------------------------------------------------- /* #

function set_iex!(b::B2_type)::B2_type

    b.temp.Ie = zeros(Int64,0)
    for p in 1:nworkers()
        ie = findmax(b.data.mp.π[b.data.kidx[p]].*(b.temp.θu[b.data.kidx[p]].-b.temp.θl[b.data.kidx[p]]))[2]
        push!(b.temp.Ie,b.data.iidx[p][ie])
    end
    
    return b
end

# */ --------------------------------------------------------------------------------------------- /* #
# */ --- lower level functions : everywhere ------------------------------------------------------ /* #
# */ --------------------------------------------------------------------------------------------- /* #

@everywhere function solv0_exact(i::Int64)::Tuple{Float64,Array{Float64,1}}

    S.ex.objf = S.ex.vars'*vcat([1.],S.temp.c[:,i])
    @objective(S.ex.m, Min, S.ex.objf)
    fix.(S.ex.m[:x],S.temp.x[:,i];force=true)
    optimize!(S.ex.m)
    θ = objective_value(S.ex.m)
    λ = dual.(FixRef.(S.ex.m[:x]))

    return θ,λ
end

@everywhere function solv1_exact(i::Int64)::Tuple{Float64,Array{Float64,1},Array{Float64,1}}

    S.ex.objf = S.ex.vars'*vcat([1.],S.temp.c[:,i])
    @objective(S.ex.m, Min, S.ex.objf)
    fix.(S.ex.m[:x],S.temp.x[:,i];force=true)
    optimize!(S.ex.m)
    θl = objective_value(S.ex.m)
    λl = dual.(FixRef.(S.ex.m[:x]))
    ϕu = value.(S.ex.m[:ϕ])

    return θl,λl,ϕu
end

@everywhere function update_m_!(s::S2_type)::S2_type

    s.m.θ[s.m.n+1:s.m.n+m0.n]   .= m0.θ[1:m0.n]
    s.m.λ[s.m.n+1:s.m.n+m0.n,:] .= m0.λ[1:m0.n,:]
    s.m.ϕ[s.m.n+1:s.m.n+m0.n,:] .= m0.ϕ[1:m0.n,:]
    s.m.x[s.m.n+1:s.m.n+m0.n,:] .= m0.x[1:m0.n,:]
    s.m.c[s.m.n+1:s.m.n+m0.n,:] .= m0.c[1:m0.n,:]
    s.m.n += m0.n

    return s
end

@everywhere function update_lb!(s::S2_type)::S2_type

    for n in s.olb.n+1:s.m.n
        s.olb.cons[n] = @constraint(s.olb.m, s.olb.m[:ϕ] + s.olb.m[:γ]'*s.m.c[n,:] >= 0. )
    end
    s.olb.help[1:s.m.n] .= s.m.θ[1:s.m.n] .- sum!(ones(s.m.n),s.m.λ[1:s.m.n,:].*s.m.x[1:s.m.n,:])
    s.olb.n = s.m.n

    return s
end

@everywhere function update_ub!(s::S2_type)::S2_type

    for n in s.oub.n+1:s.m.n
        s.oub.cons[n] = @constraint(s.oub.m, s.oub.m[:ϕ] + s.oub.m[:γ]'*s.m.x[n,:] <= 0. )
    end
    s.oub.help[1:s.m.n] .= s.m.θ[1:s.m.n] .- sum!(ones(s.m.n),s.m.ϕ[1:s.m.n,:].*s.m.c[1:s.m.n,:])
    s.oub.n = s.m.n

    return s
end

@everywhere function update_s!(s::S2_type)::S2_type

    update_m_!(s)
    update_lb!(s)
    update_ub!(s)

    return s
end

@everywhere function run_lb(i::Int64)::Tuple{Float64,Array{Float64,1}}

    S.olb.objf = S.olb.vars'*vcat([1.],S.temp.c[:,i])
    @objective(S.olb.m, Min, S.olb.objf)
    set_normalized_rhs.(S.olb.cons[1:S.olb.n], S.olb.help[1:S.olb.n].+ S.m.λ[1:S.olb.n,:]*S.temp.x[:,i])
    optimize!(S.olb.m)
    θl = objective_value(S.olb.m)
    λl = S.m.λ[1:S.m.n,:]'*dual.(S.olb.cons[1:S.olb.n])

    return θl,λl
end

@everywhere function run_ub(i::Int64)::Tuple{Float64,Array{Float64,1}}

    S.oub.objf = S.oub.vars'*vcat([1.],S.temp.x[:,i])
    @objective(S.oub.m, Max, S.oub.objf)
    set_normalized_rhs.(S.oub.cons[1:S.oub.n], S.oub.help[1:S.oub.n].+ S.m.ϕ[1:S.oub.n,:]*S.temp.c[:,i])
    optimize!(S.oub.m)
    θu = objective_value(S.oub.m)
    ϕu = S.m.ϕ[1:S.m.n,:]'*dual.(S.oub.cons[1:S.oub.n])

    return θu,ϕu
end

@everywhere function run_oracles(i::Int64)::Tuple{Float64,Float64,Array{Float64,1},Array{Float64,1}}

    θl,λl = run_lb(i)
    θu,ϕu = run_ub(i)

    return θl,θu,λl,ϕu
end


# */ --------------------------------------------------------------------------------------------- /* #

function step_a!(b::B1_type)::B1_type

    optimize!(b.rmp.m)
    b.temp.β .= value.(b.rmp.m[:β])
    b.temp.x .= value.(b.rmp.m[:x])
    for w in workers()
        @spawnat(w, Core.eval(Main, Expr(:(.=), :(S.temp.x),b.temp.x)))
    end

    return b
end

function step_a!(b::B2_type)::B2_type

    optimize!(b.rmp.m)
    b.temp.β .= value.(b.rmp.m[:β])
    b.temp.x .= value.(b.rmp.m[:x])
    for w in workers()
        @spawnat(w, Core.eval(Main, Expr(:(.=), :(S.temp.x),b.temp.x)))
    end

    return b
end

function step_b!(b::B1_type)::B1_type

    b.hist.L[b.hist.k] = exp10(6)*objective_value(b.rmp.m)

    return b
end

function step_b!(b::B2_type)::B2_type

    b.hist.L[b.hist.k] = exp10(6)*objective_value(b.rmp.m)

    return b
end

function step_c!(b::B1_type)::B1_type

    b.temp.re .= pmap(solv0_exact,1:b.data.unc.ni)
    for i in 1:b.data.unc.ni
        b.temp.θ[i] = b.temp.re[i][1]      
        b.temp.λ[:,i] .= b.temp.re[i][2]      
    end
    
    return b
end

function step_c!(b::B2_type)::B2_type

    set_iex!(b)
    b.temp.re .= pmap(solv1_exact,b.temp.Ie)
    for w in 1:nworkers()
        b.m0.θ[w]    = b.temp.re[w][1]  
        b.m0.λ[w,:] .= b.temp.re[w][2]
        b.m0.ϕ[w,:] .= b.temp.re[w][3]
        b.m0.x[w,:] .= b.temp.x[:,b.temp.Ie[w]]
        b.m0.c[w,:] .= b.temp.c[:,b.temp.Ie[w]]
    end
    sendto(procs(),m0=b.m0)
    @everywhere update_s!(S)

    return b
end

function step_d!(b::B1_type)::B1_type

    return b
end

function step_d!(b::B2_type)::B2_type

    b.temp.ro .= pmap(run_oracles,1:b.data.unc.ni)
    for i in 1:b.data.unc.ni
        b.temp.θl[i]    = b.temp.ro[i][1]
        b.temp.θu[i]    = b.temp.ro[i][2]
        b.temp.λl[:,i] .= b.temp.ro[i][3]
        b.temp.ϕu[:,i] .= b.temp.ro[i][4]
    end

    return b
end

function step_e!(b::B1_type)::B1_type

    if (b.hist.k==1)
        b.hist.U[b.hist.k] = exp10(6)*(value(b.rmp.m[:f])+b.data.mp.κ*b.data.mp.π'*b.temp.θ)
    else
        b.hist.U[b.hist.k] = min(b.hist.U[b.hist.k-1],exp10(6)*(value(b.rmp.m[:f])+b.data.mp.κ*b.data.mp.π'*b.temp.θ))
    end

    return b
end


function step_e!(b::B2_type)::B2_type

    if (b.hist.k==1)
        b.hist.U[b.hist.k] = exp10(6)*(value(b.rmp.m[:f])+b.data.mp.κ*b.data.mp.π'*b.temp.θu)
    else
        b.hist.U[b.hist.k] = min(b.hist.U[b.hist.k-1],exp10(6)*(value(b.rmp.m[:f])+b.data.mp.κ*b.data.mp.π'*b.temp.θu))
    end

    return b
end

function step_f!(b::B1_type)::B1_type

    b.temp.Δ = (b.hist.U[b.hist.k]-b.hist.L[b.hist.k])/b.hist.U[b.hist.k]*100
    if b.temp.Δ > b.data.δ
        for i in 1:b.data.unc.ni
            @constraint(b.rmp.m,b.rmp.m[:β][i] >= b.temp.θ[i]+b.temp.λ[:,i]'*(b.rmp.m[:x][:,i].-b.temp.x[:,i]))
        end
    end

    return b
end

function step_f!(b::B2_type)::B2_type

    b.temp.Δ = (b.hist.U[b.hist.k]-b.hist.L[b.hist.k])/b.hist.U[b.hist.k]*100
    if b.temp.Δ > b.data.δ
        for i in 1:b.data.unc.ni
            @constraint(b.rmp.m,b.rmp.m[:β][i] >= b.temp.θl[i]+b.temp.λl[:,i]'*(b.rmp.m[:x][:,i].-b.temp.x[:,i]))
        end
    end

    return b
end

# */ --------------------------------------------------------------------------------------------- /* #

function step_0!(b::B2_type)::B2_type

    S.temp.x[:,1] .= round.(vcat(minimum(b.data.mp.xh,dims=2)[:]*.99,minimum(b.data.unc.h,dims=1)[:].*[.99,1.01]); digits=4)
    S.temp.c[:,1] .= round.(minimum(b.data.unc.c,dims=1)[:]*.99; digits=4)
    r = solv1_exact(1)
    b.m0.n = 1
    b.m0.θ[b.m0.n] = r[1]
    b.m0.λ[b.m0.n,:] .= r[2]
    b.m0.ϕ[b.m0.n,:] .= r[3]
    b.m0.x[b.m0.n,:] .= S.temp.x[:,1]
    b.m0.c[b.m0.n,:] .= S.temp.c[:,1]
    sendto(procs(),m0=b.m0)
    @everywhere update_s!(S)
    for w in workers()
        @spawnat(w, Core.eval(Main, Expr(:(.=), :(S.temp.x),vcat(b.data.mp.xh,b.data.unc.h'))))
        @spawnat(w, Core.eval(Main, Expr(:(.=), :(S.temp.c),b.data.unc.c')))
    end
    b.m0.n = nworkers()
    b.temp.ro .= pmap(run_oracles,1:b.data.unc.ni)
    for i in 1:b.data.unc.ni
        b.temp.θl[i]    = b.temp.ro[i][1]
        b.temp.θu[i]    = b.temp.ro[i][2]
        b.temp.λl[:,i] .= b.temp.ro[i][3]
        b.temp.ϕu[:,i] .= b.temp.ro[i][4]
    end

    return b
end

function step!(b::B1_type)::B1_type

    b.hist.k += 1
    b.hist.T[b.hist.k,1] = @elapsed step_a!(b)
    b.hist.T[b.hist.k,2] = @elapsed step_b!(b)
    b.hist.T[b.hist.k,3] = @elapsed step_c!(b)
    b.hist.T[b.hist.k,4] = @elapsed step_d!(b)
    b.hist.T[b.hist.k,5] = @elapsed step_e!(b)
    b.hist.T[b.hist.k,6] = @elapsed step_f!(b)

    return b
end

function step!(b::B2_type)::B2_type

    b.hist.k += 1
    b.hist.T[b.hist.k,1] = @elapsed step_a!(b)
    b.hist.T[b.hist.k,2] = @elapsed step_b!(b)
    b.hist.T[b.hist.k,3] = @elapsed step_c!(b)
    b.hist.T[b.hist.k,4] = @elapsed step_d!(b)
    b.hist.T[b.hist.k,5] = @elapsed step_e!(b)
    b.hist.T[b.hist.k,6] = @elapsed step_f!(b)

    return b
end

# */ --------------------------------------------------------------------------------------------- /* #

function solve_Benders!(b::B1_type;print_output::Int64=1)::B1_type
    
    print_output == 1 ? print0(b) : nothing
    while (b.hist.k < b.data.J)
        step!(b)
        print_output == 1 ? print_info(b) : nothing
        (b.temp.Δ <= b.data.δ) ? break : nothing
    end
    print_output == 1 ? print_summary(b) : nothing
    
    return b
end

function solve_Benders!(b::B2_type;print_output::Int64=1)::B2_type
    
    print_output == 1 ? print0(b) : nothing
    step_0!(b)
    while (b.hist.k < b.data.J)
        step!(b)
        print_output == 1 ? print_info(b) : nothing
        (b.temp.Δ <= b.data.δ) ? break : nothing
    end
    print_output == 1 ? print_summary(b) : nothing
    
    return b
end

# */ --------------------------------------------------------------------------------------------- /* #