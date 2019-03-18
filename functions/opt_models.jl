function circ(h::Int64,H::UnitRange{Int64})::Int64
    (h==minimum(H)) ? idx = maximum(H) : idx = h-1
    return idx
end

function RMP_model(Ms::Ms_type,
                   Mp::Mp_type,
                    U::U_type )::JuMP.Model

    m = Model(solver=GurobiSolver(OutputFlag=0))

    # problem structure
    # min   f + ∑ᵢ πᵢ βᵢ)
    # s.t.  f = q₀ᵀx₀ + ∑ᵢ qᵢᵀxᵢ
    #      (x₀,x₁,..,xᵢ,..,xₙ) ∈ 𝕏
    #       βᵢ ∈ Θᵢ,  ∀i

    ##### VARIABLES ################################
    # variables of type -> f
    @variable(m, f )
    # variables of type -> x₀
    @variable(m, x0[Ms.P,Ms.I0] >= .0 )
    # variables of type -> xᵢ
    @variable(m, x[1:U.nx,Ms.I])
    # variables of type -> βᵢ
    @variable(m, β[Ms.I] >= .0 )
    ################################################

    ##### OBJ FUNCTION #############################
    # objective function -> f + ∑ᵢ πᵢβᵢ
    @objective(m, :Min, f + Mp.κ*sum(Mp.πi[i]*β[i] for i in Ms.I) )
    ################################################

    ##### CONSTRAINTS ##############################
    # constraints of type -> cx = q₀ᵀx₀ + ∑ᵢ qᵢᵀxᵢ
    @constraint(m, f_con, f >= sum(Mp.πi0[i0]*sum(Mp.c_inv[p,i0]*x0[p,i0] for p in Ms.P) for i0 in Ms.I0) + Mp.κ*sum(Mp.πi[i]*sum(Mp.c_fix[p]*x[p,i] for p in Ms.P) for i in Ms.I) )
    # constraints of type -> (x₀,x₁,..,xᵢ,..,xₙ) ∈ 𝕏
    @constraint(m, acc_x[p=Ms.P,i=Ms.I], x[p,i] == Mp.x_hist[p,i] + sum(x0[p,i0] for i0 in Mp.itoi0[i]) )
    @constraint(m, x_min[p=Ms.P,i=Ms.I], x[p,i] >= .0 )
    @constraint(m, x_max[p=Ms.P,i=Ms.I], x[p,i] <= Mp.x_max[p] )
    @constraint(m, x_νE[        i=Ms.I], x[U.nx0+1,i] ==  U.H[i,1] )
    @constraint(m, x_νD[        i=Ms.I], x[U.nx0+2,i] == -U.H[i,2] )
    # constraints of type -> βᵢ ∈ Θᵢ,  ∀i
    ################################################

    return m
end

function SP_model(Ps::Ps_type,
                  Pp::Pp_type,
                   U::U_type )::JuMP.Model

    c,xf = zeros(U.nc),zeros(U.nx)
    ξ = vcat(ones(U.nx0).*0.001,ones(2))

    m = Model(solver=GurobiSolver(OutputFlag=0,Method=2))

    # problem structure
    # min   c₀ + cᵀϕ
    # s.t.  c₀  = C₀y
    #       ϕ   = C y
    #       y ∈ 𝕐
    #       A y ≦ B x
    #       x = xf :(λ)

    ##### VARIABLES ################################
    # variables of type -> ϕ
    @variable(m, ϕ[1:U.nc] ) # (MWh)
    # variables of type -> c₀
    @variable(m, c0 ) # (£)
    # variables of type -> y
    @variable(m, yG[Ps.G,Ps.S,Ps.H] >= 0  ) # power generation from conventional unit g at hour h of season s (MW)
    @variable(m, yI[Ps.B,Ps.S,Ps.H] >= 0  ) # charging power of storage unit b at hour h of season s (MW)
    @variable(m, yO[Ps.B,Ps.S,Ps.H] >= 0  ) # discharging power of storage unit b at hour h of season s (MW)
    @variable(m, yL[Ps.B,Ps.S,Ps.H] >= 0  ) # energy level of storage unit b at hour h of season s (MWh)
    @variable(m, yS[     Ps.S,Ps.H] >= 0  ) # shedded demand at hour h of season s (MW)
    # variables of type -> x
    @variable(m, x[1:U.nx] ) # installed capacity of technology p (MW)
    ################################################

    ##### OBJ FUNCTION #############################
    # objective function -> c₀ + cᵀϕ
    @objective(m, :Min, c0 + sum(c[n]*ϕ[n] for n in 1:U.nc))
    ################################################

    ##### CONSTRAINTS ##############################
    # constraints of type -> c₀  = C₀y
    @constraint(m, c0_con, c0 == sum(sum(sum((Pp.c_OMvarG[g]+Pp.c_fuelG[g]/Pp.ηG[g])*yG[g,s,h] for g in 1:5) + Pp.c_OMvarG[6]*yG[6,s,h] + Pp.c_shed*yS[s,h] for h in Ps.H)*Pp.α for s in Ps.S))
    # constraints of type -> ϕ  = C y
    @constraint(m, ϕ1_con , ϕ[1] == sum(sum(sum(Pp.em_co2G[g]/Pp.ηG[g]*yG[g,s,h] for g in Ps.G) for h in Ps.H)*Pp.α for s in Ps.S))
    @constraint(m, ϕ2_con , ϕ[2] == sum(sum(                1/Pp.ηG[6]*yG[6,s,h]                for h in Ps.H)*Pp.α for s in Ps.S))
    # constraints of type -> y ∈ 𝕐
    @constraint(m, Y_con[b=Ps.B,s=Ps.S,h=Ps.H], yL[b,s,h] - yL[b,s,circ(h,Ps.H)] == Pp.ηB[b-Ps.b0]*yI[b,s,h] - yO[b,s,h])
    # constraints of type -> A y ≦ B x
    @constraint(m, A1_con[g=Ps.G,s=Ps.S,h=Ps.H], yG[g,s,h] - yG[g,s,circ(h,Ps.H)] <= Pp.rampG[g]*x[g]) # impose ramp up constraint of generator g
    @constraint(m, A2_con[g=Ps.G,s=Ps.S,h=Ps.H], yG[g,s,circ(h,Ps.H)] - yG[g,s,h] <= Pp.rampG[g]*x[g]) # impose ramp dw constraint of generator g
    @constraint(m, A3_con[g=Ps.G,s=Ps.S,h=Ps.H], yG[g,s,h] <= x[g]) # impose capacity of generator g
    @constraint(m, A4_con[b=Ps.B,s=Ps.S,h=Ps.H], yI[b,s,h] <= Pp.PB[b-Ps.b0]*x[b]) # impose charging capacity of storage b
    @constraint(m, A5_con[b=Ps.B,s=Ps.S,h=Ps.H], yO[b,s,h] <= Pp.PB[b-Ps.b0]*x[b]) # impose discharging capacity of storage b
    @constraint(m, A6_con[b=Ps.B,s=Ps.S,h=Ps.H], yL[b,s,h] <= x[b]) # impose capacity of storage b
    @constraint(m, A7_con, sum(yG[g,s,h]*Pp.em_co2G[g]/Pp.ηG[g] for g in Ps.G for s in Ps.S for h = Ps.H) <= Pp.co2_lim*x[U.nx0+1]) # impose emission limit
    @constraint(m, A8_con[s=Ps.S,h=Ps.H], sum(yG[g,s,h] for g in Ps.G) + sum(yO[b,s,h]-yI[b,s,h] for b in Ps.B) + yS[s,h] >= -x[U.nx0+2]*Pp.PD[s,h] - sum(Pp.PR[r-Ps.r0,s,h]*x[r] for r in Ps.R)) # impose energy balance
    # constraints of type -> x = xf :(λ)
    @constraint(m, λ[n=1:U.nx], ξ[n]*x[n] == xf[n])
    ################################################

    return m
end

function Oracle_LB(U::U_type)::JuMP.Model

         m = Model(solver=GurobiSolver(OutputFlag=0))

         # problem structure
         # min   ϕ + γᵀcᵢ
         # s.t.  ϕ + γᵀcₛ ≧ θₛ + λᵀₛ(x-xₛ),   ∀s
         #       x = xᵢ

         ##### VARIABLES ################################
         # variables of type -> ϕ
         @variable(m, ϕ)
         # variables of type -> γ
         @variable(m, γ[1:U.nc] >= .0)
         # variables of type -> x
         @variable(m, x[1:U.nx])
         ################################################

         ##### OBJ FUNCTION #############################
         # objective function -> ϕ + γᵀcᵢ
         @objective(m, :Min, 0. )
         ################################################

         ##### CONSTRAINTS ##############################
         # constraints of type -> ϕ + γᵀcₛ ≧ θₛ + λᵀₛ(x-xₛ),   ∀s
         # constraints of type -> x = xᵢ
         @constraint(m,x_fix[n=1:U.nx], x[n] == 0.)
         ################################################
         return m
end

function Oracle_UB(U::U_type)::JuMP.Model

         m = Model(solver=GurobiSolver(OutputFlag=0))

         # problem structure
         # max   ϕ + γᵀxᵢ
         # s.t.  ϕ + γᵀxₛ ≦ θₛ + ϕᵀₛ(c-cₛ),   ∀s
         #       c = cᵢ

         ##### VARIABLES ################################
         # variables of type -> ϕ
         @variable(m, ϕ)
         # variables of type -> γ
         @variable(m, γ[1:U.nx] <= 0.)
         # variables of type -> c
         @variable(m, c[1:U.nc])
         ################################################

         ##### OBJ FUNCTION #############################
         # objective function -> ϕ + γᵀxᵢ
         @objective(m, :Max, 0.)
         ################################################

         ##### CONSTRAINTS ##############################
         # constraints of type -> ϕ + γᵀxₛ ≦ θₛ + ϕᵀₛ(c-cₛ),   ∀s
         # constraints of type -> x = xᵢ
         @constraint(m,c_fix[n=1:U.nc], c[n] == 0.)
         ################################################
         return m
end
