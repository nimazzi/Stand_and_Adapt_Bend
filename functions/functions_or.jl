mutable struct J_type_or
        x::Array{Float64,2}
        β::Vector{Float64}
   θunder::Vector{Float64}
   λunder::Array{Float64,2}
    θover::Vector{Float64}
        Δ::Float64
end

mutable struct B_type_or
    L::Vector{Float64}
    U::Vector{Float64}
    t::Vector{Float64}
   ta::Vector{Float64}
   tb::Vector{Float64}
   tc::Vector{Float64}
   td::Vector{Float64}
end

mutable struct S_type_or
    θ::Vector{Float64}
    λ::Array{Float64,2}
    ϕ::Array{Float64,2}
    x::Array{Float64,2}
    c::Array{Float64,2}
end

mutable struct R_type_or
     m::JuMP.Model
  vars::Array{JuMP.Variable,2}
  cons::Array{JuMP.ConstraintRef,2}
end

mutable struct E_type_or
     m::JuMP.Model
  vars::Array{JuMP.Variable,1}
end

mutable struct LB_oracle
         m::JuMP.Model
      cons::Array{JuMP.ConstraintRef,1}
      vars::Array{JuMP.Variable,1}
end

mutable struct UB_oracle
           m::JuMP.Model
        cons::Array{JuMP.ConstraintRef,1}
        vars::Array{JuMP.Variable,1}
end

mutable struct O_type_or
          LB::LB_oracle
          UB::UB_oracle
end

mutable struct T_type_or
       LBcon::Array{Float64,1}
       UBcon::Array{Float64,1}
end


function generate_structs_or(Ms::Ms_type,
                             Mp::Mp_type,
                             Ps::Ps_type,
                             Pp::Pp_type,
                              U::U_type)::Tuple{R_type_or,E_type_or,O_type_or,J_type_or,B_type_or,S_type_or,T_type_or}

     m = RMP_model(Ms,Mp,U)
     vars = Array{JuMP.Variable,2}(undef,U.nx+1,U.nI)
     for i in 1:U.nI vars[:,i] = vcat([m[:β][i]],m[:x][:,i]) end
     @constraintref cons[1:ITmax,1:U.nI]
     #   R_type_or(m,cons)
     R = R_type_or(m,vars,cons)
     solve(R.m)

     m = SP_model(Ps,Pp,U)
     vars = vcat(m[:c0],m[:ϕ][:])
     #   E_type_or(m,vars)
     E = E_type_or(m,vars)
     solve(E.m)

     m = Oracle_LB(U); solve(m)
     @constraintref cons[1:ITmax]
     vars = vcat([m[:ϕ]],[m[:γ][n] for n in 1:U.nc])
     #    LB_oracle(m,cons,c_vars,o_vars)
     LB = LB_oracle(m,cons,vars)

     m = Oracle_UB(U); solve(m)
     @constraintref cons[1:ITmax]
     vars = vcat([m[:ϕ]],[m[:γ][n] for n in 1:U.nx])
     #    UB_oracle(m,cons,c_vars,o_vars)
     UB = UB_oracle(m,cons,vars)

     #   O_type_or(LB,UB)
     O = O_type_or(LB,UB)

     #   J_type_or(x                                      ,β          ,θunder     ,λunder          ,θover      ,Δ )
     J = J_type_or(round.(getvalue(R.m[:x][:,:]);digits=4),zeros(U.nI),zeros(U.nI),zeros(U.nx,U.nI),zeros(U.nI),0.)

     #   B_type_or(L       ,U       ,t       ,ta      ,tb      ,tc      ,td      )
     B = B_type_or(zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0))

     #   S_type_or(θ       ,λ            ,ϕ            ,x            ,c            )
     S = S_type_or(zeros(0),zeros(0,U.nx),zeros(0,U.nc),zeros(0,U.nx),zeros(0,U.nc))

     #   T_type_or(LBcon   ,UBcon   )
     T = T_type_or(zeros(1),zeros(1))

    return R,E,O,J,B,S,T
end

function do_step0_or(Ms::Ms_type,
                     Mp::Mp_type,
                      U::U_type,
                      E::E_type_or,
                      S::S_type_or,
                      O::O_type_or,
                      J::J_type_or,
                      T::T_type_or)::Tuple{E_type_or,S_type_or,O_type_or,J_type_or,T_type_or}

     x_min = vcat([minimum(Mp.x_hist[p,:])*.99 for p in Ms.P],[minimum(U.H[:,1])*.99],-maximum(U.H[:,2])*1.01); x_min = round.(x_min;digits=4)
     c_min = vcat([minimum(U.C[:,1])*.99],[minimum(U.C[:,2])*.99]); c_min = round.(c_min;digits=4)

     @objective(E.m, :Min, AffExpr(E.vars,vcat([1.],c_min),.0))
     JuMP.setRHS.(E.m[:λ][:],x_min[:])
     solve(E.m)

     S.θ = vcat(S.θ,getobjectivevalue(E.m)/10^8)
     S.λ = vcat(S.λ,getdual.(E.m[:λ][:])'./10^8)
     S.ϕ = vcat(S.ϕ,getvalue.(E.m[:ϕ][:])'./10^8)
     S.x = vcat(S.x,x_min')
     S.c = vcat(S.c,c_min')
     ns  = length(S.θ)

     O.LB.cons[ns] = @constraint(O.LB.m, AffExpr(O.LB.vars,vcat([1.],S.c[ns,:]),0.) >= 0.)
     T.LBcon .= S.θ .- reshape(sum(S.λ.*S.x,dims=2),ns)
     for i in 1:U.nI
          # lower bound
          @objective(O.LB.m, :Min, AffExpr(O.LB.vars,vcat([1.],U.C[i,:]),.0))
          JuMP.setRHS.(O.LB.cons[1:ns],T.LBcon .+ S.λ*J.x[:,i])
          solve(O.LB.m)
          J.θunder[i] = getobjectivevalue(O.LB.m)*10^8
          J.λunder[:,i] .= S.λ'*getdual.(O.LB.cons[1:ns]).*10^8
     end

     O.UB.cons[ns] = @constraint(O.UB.m, AffExpr(O.UB.vars,vcat([1.],S.x[ns,:]),0.) <= 0.)
     T.UBcon .= S.θ .- reshape(sum(S.ϕ.*S.c,dims=2),ns)
     for i in 1:U.nI
          # upper bound
          @objective(O.UB.m, :Max, AffExpr(O.UB.vars,vcat([1.],J.x[:,i]),.0))
          JuMP.setRHS.(O.UB.cons[1:ns],T.UBcon .+ S.ϕ*U.C[i,:])
          solve(O.UB.m)
          J.θover[i] = getobjectivevalue(O.UB.m)*10^8
     end

    return E,S,O,J,T
end

function do_step_A_or(Mp::Mp_type,
                       R::R_type_or,
                       J::J_type_or,
                       B::B_type_or)::Tuple{R_type_or,J_type_or,B_type_or,Int64}

    solve(R.m)
    J.β .= getvalue(R.m[:β][:])
    J.x .= round.(getvalue(R.m[:x][:,:]);digits=4)
    push!(B.L,getobjectivevalue(R.m))
    ω = findmax(Mp.πi.*(J.θover.-J.θunder))[2]

    return R,J,B,ω
end

function do_step_B_or(U::U_type,
                      ω::Int64,
                      J::J_type_or,
                      E::E_type_or,
                      S::S_type_or)::Tuple{E_type_or,S_type_or}

    @objective(E.m, :Min, AffExpr(E.vars,vcat([1.],U.C[ω,:]),.0))
    JuMP.setRHS.(E.m[:λ][:],J.x[:,ω])
    solve(E.m)

    S.θ = vcat(S.θ,getobjectivevalue(E.m)/10^8)
    S.λ = vcat(S.λ,getdual.(E.m[:λ][:])'./10^8)
    S.ϕ = vcat(S.ϕ,getvalue.(E.m[:ϕ][:])'./10^8)
    S.x = vcat(S.x,J.x[:,ω]')
    S.c = vcat(S.c,U.C[ω,:]')

    return E,S
end

function do_step_C_or(S::S_type_or,
                      U::U_type,
                      O::O_type_or,
                      J::J_type_or,
                      T::T_type_or)::Tuple{O_type_or,J_type_or,T_type_or}

    ns=length(S.θ)

    push!(T.LBcon,0.)
    O.LB.cons[ns] = @constraint(O.LB.m, AffExpr(O.LB.vars,vcat([1.],S.c[ns,:]),0.) >= 0.)
    T.LBcon .= S.θ .- reshape(sum(S.λ.*S.x,dims=2),ns)
    for i in 1:U.nI
         # lower bound
         @objective(O.LB.m, :Min, AffExpr(O.LB.vars,vcat([1.],U.C[i,:]),.0))
         JuMP.setRHS.(O.LB.cons[1:ns],T.LBcon .+ S.λ*J.x[:,i])
         solve(O.LB.m)
         J.θunder[i] = getobjectivevalue(O.LB.m)*10^8
         J.λunder[:,i] .= S.λ'*getdual.(O.LB.cons[1:ns]).*10^8
    end

    push!(T.UBcon,0.)
    O.UB.cons[ns] = @constraint(O.UB.m, AffExpr(O.UB.vars,vcat([1.],S.x[ns,:]),0.) <= 0.)
    T.UBcon .= S.θ .- reshape(sum(S.ϕ.*S.c,dims=2),ns)
    for i in 1:U.nI
         # upper bound
         @objective(O.UB.m, :Max, AffExpr(O.UB.vars,vcat([1.],J.x[:,i]),.0))
         JuMP.setRHS.(O.UB.cons[1:ns],T.UBcon .+ S.ϕ*U.C[i,:])
         solve(O.UB.m)
         J.θover[i] = getobjectivevalue(O.UB.m)*10^8
    end

    return O,J,T
end

function do_step_D_or(Mp::Mp_type,
                       J::J_type_or,
                       U::U_type,
                       R::R_type_or,
                       B::B_type_or)::Tuple{R_type_or,B_type_or}

    (length(B.U)==0) ? push!(B.U,getvalue(R.m[:f])+Mp.κ*sum(Mp.πi.*J.θover)) : push!(B.U,min(B.U[end],getvalue(R.m[:f])+Mp.κ*sum(Mp.πi.*J.θover)))
    J.Δ = (B.U[end]-B.L[end])/B.U[end]*100.
    for i in 1:U.nI R.cons[length(B.L),i] = @constraint(R.m, AffExpr(R.vars[:,i],vcat([1.],-J.λunder[:,i]),.0) >= J.θunder[i]-sum(J.λunder[:,i].*J.x[:,i])) end

    return R,B
end

function B_time_or(B::B_type_or,
                  ta::Float64,
                  tb::Float64,
                  tc::Float64,
                  td::Float64)::B_type_or

    push!(B.ta,ta         )
    push!(B.tb,   tb)
    push!(B.tc,      tc)
    push!(B.td,         td)
    push!(B.t ,ta+tb+tc+td)

    return B
end

function do_step_or(Mp::Mp_type,
                     U::U_type,
                     R::R_type_or,
                     E::E_type_or,
                     O::O_type_or,
                     S::S_type_or,
                     B::B_type_or,
                     J::J_type_or,
                     T::T_type_or)::Tuple{R_type_or,E_type_or,O_type_or,S_type_or,B_type_or,J_type_or,T_type_or}

    ta = @elapsed R,J,B,ω = do_step_A_or(Mp,R,J,B)
    tb = @elapsed E,S     = do_step_B_or(U,ω,J,E,S)
    tc = @elapsed O,J,T   = do_step_C_or(S,U,O,J,T)
    td = @elapsed R,B     = do_step_D_or(Mp,J,U,R,B)

    B = B_time_or(B,ta,tb,tc,td)

    return R,E,O,S,B,J,T
end

function print_init_or(case::Int64,
                         Ms::Ms_type)

    println(" ")
    println("decomposition algorithm:")
    println("*" ^ 50)
    println(" algorithm          : " * "Benders with adaptive oracles")
    println(" case               : $(case)" )
    println(" investment  nodes  : $(maximum(Ms.I0))"   )
    println(" operational nodes  : $(maximum(Ms.I ))"    )
    println("-" ^ 50)

end

function string_fcn_ABC_or(j::Int64,
                           Δ::Float64,
                           t::Float64)::String

    str_A = " j = " * (" " ^ (4-length("$(j)"))) * "$(j), "
    B = "$(round(Δ;digits=3))" * "0" ^ (5 - length("$(round(Δ%1;digits=3))"))
    str_B = "Δ = " * (" " ^ (8-length(B))) * "$(B) %, "
    C = "$(round(t;digits=1))"
    str_C = "t = " * (" " ^ (7-length(C))) * "$(C) s"

    str_ABC = str_A * str_B * str_C

    return str_ABC
end

function print_info_or(J::J_type_or,
                       B::B_type_or)

    j = length(B.L)
    Δ = J.Δ
    t = B.t[j]
    str = string_fcn_ABC_or(j,Δ,t)
    println(str)

end

function str_fcn1_or(x_ay::Vector{Int64},
                     tech::String)::String

    I = length(x_ay)
    str = " " * tech * " " ^ (9-length(tech))
    for i in 1:I
        x = string(x_ay[i])
        str *= (" " ^ (7-length(x)) * x)
    end

    return str
end

function print_end_summary_or(Ms::Ms_type,
                               U::U_type,
                               R::R_type_or,
                               B::B_type_or,
                            case::Int64)

    println(""); println("*" ^ 75); println("*" ^ 75)
    x0 = zeros(U.nx0,maximum(Ms.I0)); x0 .= getvalue.(R.m[:x0][:,:])
    x1_lo,x2_lo = convert(Array{Int64,1},round.(x0[:,1];digits=0)),convert(Array{Int64,2},round.(x0[:,2:end];digits=0))
    U_lo,L_lo,t_lo,tA_lo,tB_lo,tC_lo = B.U,B.L,B.t,B.ta,B.tb,B.tc
    str_tech = ["coal","coalccs","OCGT","CCGT","diesel","nuclear","pumpL","pumpH","lithium","onwind","offwind","solar"]
    int_A = " tech.    "
    for i in 1:length(x1_lo[1,:])
        st  = "ω$i"
        int_A *= " " ^ (7 - length(st)) * st
    end
    if (case < 4)
        int_B = " tech.    "
        for i in 1:length(x2_lo[1,:]) st  = "ω$i"; int_B *= " " ^ (7 - length(st)) * st; end
        Ni,No = 3^(case -1)+1,3^(2case -2)+3^(case -1)
        obj = round(U_lo[end] / 1e11;digits=3)
        ebin,cbin,ubin = 0,0,0
        if (case >= 2) ebin = 1 end; if (case >= 3) cbin = 1 end; if (case >= 4) ubin = 1 end
        (ebin == 0) ? estr = "deterministic" : estr = "uncertain"
        (cbin == 0) ? cstr = "deterministic" : cstr = "uncertain"
        (ubin == 0) ? ustr = "deterministic" : ustr = "uncertain"
        println(" ")
        println(" co2 emission limit : $(estr)"); println(" co2 emission cost  : $(cstr)"); println(" uranium cost       : $(ustr)")
        println(" investment nodes   : $(Ni)")  ; println(" operational nodes  : $(No)")  ; println(" optimal objective  : $(obj) x 10^11 £")
        println(" "); println(" optimal investments @ 0 years"); println(" " * "-" ^ (length(int_A)-1))
        println(int_A); println("-" ^ (length(int_A)-1))
        for p in 1:12 println(str_fcn1_or(x1_lo[p,:],str_tech[p])) end
        println(" " * "-" ^ (length(int_A)-1)); println(" ")
        println(" optimal investments @ 5 years"); println(" " * "-" ^ (length(int_B)-1)); println(int_B); println(" " * "-" ^ (length(int_B)-1))
        for p in 1:12 println(str_fcn1_or(x2_lo[p,:],str_tech[p])) end
        println(" " * "-" ^ (length(int_B)-1)); println(" ")
    end
    if case==4
        int_B = [" tech.    "," tech.    "," tech.    "]
        for j in 1:3 for i in 1:Int64(length(x2_lo[1,:])/3) st  = "ω$(9*(j-1)+i)"; int_B[j] *= " " ^ (7 - length(st)) * st; end end
        Ni,No = 3^(case-1)+1,3^(2case-2)+3^(case-1)
        obj = round(U_lo[end] / 1e11; digits=3)
        ebin,cbin,ubin = 0,0,0
        if (case >= 2) ebin = 1 end; if (case >= 3) cbin = 1 end; if (case >= 4) ubin = 1 end
        (ebin == 0) ? estr = "deterministic" : estr = "uncertain"
        (cbin == 0) ? cstr = "deterministic" : cstr = "uncertain"
        (ubin == 0) ? ustr = "deterministic" : ustr = "uncertain"
        println(" ")
        println(" co2 emission limit : $(estr)"); println(" co2 emission cost  : $(cstr)"); println(" uranium cost       : $(ustr)")
        println(" investment nodes   : $(Ni)")  ; println(" operational nodes  : $(No)")  ; println(" optimal objective  : $(obj) x 10^11 £")
        println(" "); println(" optimal investments @ 0 years"); println(" " * "-" ^ (length(int_A)-1))
        println(int_A); println("-" ^ (length(int_A)-1))
        for p in 1:12 println(str_fcn1_or(x1_lo[p,:],str_tech[p])) end
        println(" " * "-" ^ (length(int_A)-1)); println(" ")
        println(" optimal investments @ 5 years")
        for j in 1:3
            println(" " * "-" ^ (length(int_B[j])-1)); println(int_B[j]); println(" " * "-" ^ (length(int_B[j])-1))
            for p in 1:12 println(str_fcn1_or(x2_lo[p,9*(j-1)+1:9*j],str_tech[p])) end
            println(" " * "-" ^ (length(int_B[j])-1)); println(" ")
        end
    end
    println(" computational results:")
    println(" " * "-" ^ 67)
    println(" ϵ-target (%)" * " | " * "iters     time (s)" * "   " * "RMP (%)   SP (%)   Oracles (%)")
    println(" " * "-" ^ 67)
    for ϵ in [1.00,0.10,0.01]
       Δ_lo    = (U_lo.-L_lo)./U_lo*100
       N_lo    = findmin(max.(Δ_lo .- ϵ,0))[2]
       tt_lo,ta_lo,tb_lo,tc_lo = sum(t_lo[1:N_lo]),sum(tA_lo[1:N_lo]),sum(tB_lo[1:N_lo]),sum(tC_lo[1:N_lo])
       RMPlo,SPlo,Olo = round(ta_lo/tt_lo*100;digits=2),round(tb_lo/tt_lo*100;digits=2),round(tc_lo/tt_lo*100;digits=2)
       str_ϵ   = " $(round(ϵ;digits=2))"; str_ϵ   *= "0" ^ (5-length(str_ϵ))
       str_Nlo = "$(N_lo)"  ; str_Nlo *= " " ^ (10-length(str_Nlo))
       str_tlo = "$(Int64(round(tt_lo;digits=0)))" ; str_tlo *= " " ^ (9-length(str_tlo))
       str_Rlo     = "$(RMPlo)"; str_Rlo = " " ^ (3-findfirst(isequal('.'),str_Rlo)) * str_Rlo; str_Rlo *= "0" ^ (5-length(str_Rlo)); str_Rlo *= " " ^ (10-length(str_Rlo))
       str_Slo     = "$(SPlo)" ; str_Slo = " " ^ (4-findfirst(isequal('.'),str_Slo)) * str_Slo; str_Slo *= "0" ^ (6-length(str_Slo)); str_Slo *= " " ^ (9-length(str_Slo))
       str_Olo     = "$(Olo)"  ; str_Olo = " " ^ (3-findfirst(isequal('.'),str_Olo)) * str_Olo; str_Olo *= "0" ^ (5-length(str_Olo))

       println(str_ϵ * "         | " * str_Nlo * str_tlo * "  " * str_Rlo * str_Slo * str_Olo)
    end
    println(" " * "-" ^ 67)
    println(""); println("*" ^ 75); println("*" ^ 75)
end
