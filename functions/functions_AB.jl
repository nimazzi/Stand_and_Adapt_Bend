################################################################################
###### Algorithm 2: 𝖠𝖽𝖺𝗉𝗍_𝖡𝖾𝗇𝖽 (Benders with Adaptive Oracle) #################
###### iteration (j)
# solve 𝐑𝐌𝐏 and obtain 𝐱ʲ and βʲ                       \* step A  solve master problem            *\
# set Lʲ:= 𝑓(𝐱ʲ)+∑ᵢβʲᵢ                                   \* step B  compute lower bound            *\
# set ξ:=0 and 𝓔:=𝓘                                     \* step C  solve subproblem
# repeat
#     select iex ∈ 𝓔
#     set 𝓔 := 𝓔 \ iex
#     solve 𝐒𝐏(iex) at (xʲᵢ,cᵢ) and obtain θʲᵢ, λʲᵢ, ϕʲᵢ
#     if (xʲᵢ,θʲᵢ,λʲᵢ) is non-redundant
#         set ξ:=1
#         set 𝓢:=𝓢 ∪ {(xʲᵢ,cᵢ,θʲᵢ,λʲᵢ, ϕʲᵢ)}
# until ξ=1 or 𝓔=∅                                                                                 *\
# for i ∈ 𝓘                                             \* step D  solve oracles
#     solve 𝓐LB(xʲᵢ,cᵢ) and obtain θlʲᵢ and λlʲᵢ
#     solve 𝓐UB(xʲᵢ,cᵢ) and obtain θuʲᵢ                                                           *\
# set Uʲ:=𝑓(𝐱ʲ)+∑ᵢθuʲᵢ                                   \* step E  compute upper bound            *\
# for i ∈ 𝓘                                             \* step F  update master problem
#     add cut (xʲᵢ,θlʲᵢ,λlʲᵢ) to 𝐑𝐌𝐏                                                             *\
################################################################################

mutable struct J_type_AB
        x::Array{Float64,2}
        β::Vector{Float64}
   θunder::Vector{Float64}
   λunder::Array{Float64,2}
    θover::Vector{Float64}
        Δ::Float64
end

mutable struct B_type_AB
    L::Vector{Float64}
    U::Vector{Float64}
    t::Vector{Float64}
   tA::Vector{Float64}
   tB::Vector{Float64}
   tC::Vector{Float64}
   tD::Vector{Float64}
   tE::Vector{Float64}
   tF::Vector{Float64}
end

mutable struct S_type_AB
    θ::Vector{Float64}
    λ::Array{Float64,2}
    ϕ::Array{Float64,2}
    x::Array{Float64,2}
    c::Array{Float64,2}
end

mutable struct R_type_AB
     m::JuMP.Model
  vars::Array{JuMP.Variable,2}
  cons::Array{JuMP.ConstraintRef,2}
    θℓ::Array{Array{Float64,1},1}
    xℓ::Array{Array{Float64,2},1}
    λℓ::Array{Array{Float64,2},1}
end

mutable struct E_type_AB
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

mutable struct O_type_AB
          LB::LB_oracle
          UB::UB_oracle
end

mutable struct T_type_AB
       LBcon::Array{Float64,1}
       UBcon::Array{Float64,1}
end


function gen_structs_AB(Ms::Ms_type,
                        Mp::Mp_type,
                        Ps::Ps_type,
                        Pp::Pp_type,
                         U::U_type)::Tuple{R_type_AB,E_type_AB,O_type_AB,J_type_AB,B_type_AB,S_type_AB,T_type_AB}

     m = RMP_model(Ms,Mp,U)
     vars = Array{JuMP.Variable,2}(undef,U.nx+1,U.nI)
     for i in 1:U.nI vars[:,i] = vcat([m[:β][i]],m[:x][:,i]) end
     @constraintref cons[1:ITmax,1:U.nI]
     #   R_type_AB(m,vars,cons)
     R = R_type_AB(m,vars,cons,[zeros(0) for i in 1:U.nI],[zeros(0,U.nx) for i in 1:U.nI],[zeros(0,U.nx) for i in 1:U.nI])
     solve(R.m)

     m = SP_model(Ps,Pp,U)
     vars = vcat(m[:c0],m[:ϕ][:])
     #   E_type_AB(m,vars)
     E = E_type_AB(m,vars)
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

     #   O_type_AB(LB,UB)
     O = O_type_AB(LB,UB)

     #   J_type_AB(x                                      ,β          ,θunder     ,λunder          ,θover      ,Δ )
     J = J_type_AB(round.(getvalue(R.m[:x][:,:]);digits=4),zeros(U.nI),zeros(U.nI),zeros(U.nx,U.nI),zeros(U.nI),0.)

     #   B_type_AB(L       ,U       ,t       ,tA      ,tB      ,tC      ,tD      ,tE      ,tF      )
     B = B_type_AB(zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0))

     #   S_type_AB(θ       ,λ            ,ϕ            ,x            ,c            )
     S = S_type_AB(zeros(0),zeros(0,U.nx),zeros(0,U.nc),zeros(0,U.nx),zeros(0,U.nc))

     #   T_type_AB(LBcon   ,UBcon   )
     T = T_type_AB(zeros(1),zeros(1))

    return R,E,O,J,B,S,T
end

function Adapt_Bend_step_0(Ms::Ms_type,
                           Mp::Mp_type,
                            U::U_type,
                            E::E_type_AB,
                            S::S_type_AB,
                            O::O_type_AB,
                            J::J_type_AB,
                            T::T_type_AB)::Tuple{E_type_AB,S_type_AB,O_type_AB,J_type_AB,T_type_AB}

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

function Adapt_Bend_step_A(R::R_type_AB,
                           J::J_type_AB)::Tuple{R_type_AB,J_type_AB}

    solve(R.m)
    J.β .= getvalue(R.m[:β][:])
    J.x .= round.(getvalue(R.m[:x][:,:]);digits=4)

    return R,J
end

function Adapt_Bend_step_B(R::R_type_AB,
                           B::B_type_AB)::B_type_AB

    push!(B.L,getobjectivevalue(R.m))

    return B
end

function Adapt_Bend_step_C(Mp::Mp_type,
                            U::U_type,
                           R::R_type_AB,
                           J::J_type_AB,
                           E::E_type_AB,
                           S::S_type_AB)::Tuple{E_type_AB,S_type_AB}

    ξ,ℰ = 0,[i for i in 1:U.nI]
    while true
        iex = findmax([Mp.πi[i]*(J.θover[i]-J.θunder[i]) for i in ℰ])[2]
        filter!(x->x≠iex,ℰ)
        @objective(E.m, :Min, AffExpr(E.vars,vcat([1.],U.C[iex,:]),.0))
        JuMP.setRHS.(E.m[:λ][:],J.x[:,iex])
        solve(E.m)
        θ = getobjectivevalue(E.m)
        θx = R.θℓ[iex] .+ R.λℓ[iex]*(repeat(J.x[:,iex]',outer=length(R.θℓ[iex])).-R.xℓ[iex])'
        if (any(x->x==θ,θx)==false)
            ξ=1
            S.θ = vcat(S.θ,getobjectivevalue(E.m)/10^8)
            S.λ = vcat(S.λ,getdual.(E.m[:λ][:])'./10^8)
            S.ϕ = vcat(S.ϕ,getvalue.(E.m[:ϕ][:])'./10^8)
            S.x = vcat(S.x,J.x[:,iex]')
            S.c = vcat(S.c,U.C[iex,:]')
        end
        if (ξ==1 || length(ℰ)==0)
            break
        end
    end
    return E,S
end

function Adapt_Bend_step_D(S::S_type_AB,
                           U::U_type,
                           O::O_type_AB,
                           J::J_type_AB,
                           T::T_type_AB)::Tuple{O_type_AB,J_type_AB,T_type_AB}

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

function Adapt_Bend_step_E(Mp::Mp_type,
                            J::J_type_AB,
                            R::R_type_AB,
                            B::B_type_AB)::B_type_AB

    (length(B.U)==0) ? push!(B.U,getvalue(R.m[:f])+Mp.κ*sum(Mp.πi.*J.θover)) : push!(B.U,min(B.U[end],getvalue(R.m[:f])+Mp.κ*sum(Mp.πi.*J.θover)))

    return B
end

function Adapt_Bend_step_F(J::J_type_AB,
                           U::U_type,
                           R::R_type_AB,
                           B::B_type_AB)::Tuple{R_type_AB,J_type_AB}

    for i in 1:U.nI R.cons[length(B.L),i] = @constraint(R.m, AffExpr(R.vars[:,i],vcat([1.],-J.λunder[:,i]),.0) >= J.θunder[i]-sum(J.λunder[:,i].*J.x[:,i])) end
    for i in 1:U.nI
        push!(R.θℓ[i],J.θunder[i])
        R.λℓ[i] = vcat(R.λℓ[i],J.λunder[:,i]')
        R.xℓ[i] = vcat(R.xℓ[i],J.x[:,i]')
    end
    J.Δ = (B.U[end]-B.L[end])/B.U[end]*100.

    return R,J
end

function B_time_AB(B::B_type_AB,
                  tA::Float64,
                  tB::Float64,
                  tC::Float64,
                  tD::Float64,
                  tE::Float64,
                  tF::Float64)::B_type_AB

    push!(B.tA,tA         )
    push!(B.tB,   tB)
    push!(B.tC,      tC)
    push!(B.tD,         tD)
    push!(B.tE,            tE)
    push!(B.tF,               tF)
    push!(B.t ,tA+tB+tC+tD+tE+tF)

    return B
end

function iter_AB(Mp::Mp_type,
                  U::U_type,
                  R::R_type_AB,
                  E::E_type_AB,
                  O::O_type_AB,
                  S::S_type_AB,
                  B::B_type_AB,
                  J::J_type_AB,
                  T::T_type_AB)::Tuple{R_type_AB,E_type_AB,O_type_AB,S_type_AB,B_type_AB,J_type_AB,T_type_AB}

    tA = @elapsed R,J   = Adapt_Bend_step_A(R,J)
    tB = @elapsed B     = Adapt_Bend_step_B(R,B)
    tC = @elapsed E,S   = Adapt_Bend_step_C(Mp,U,R,J,E,S)
    tD = @elapsed O,J,T = Adapt_Bend_step_D(S,U,O,J,T)
    tE = @elapsed B     = Adapt_Bend_step_E(Mp,J,R,B)
    tF = @elapsed R,J   = Adapt_Bend_step_F(J,U,R,B)

    B = B_time_AB(B,tA,tB,tC,tD,tE,tF)

    return R,E,O,S,B,J,T
end

function print_init_AB(case::Int64,
                         Ms::Ms_type)

    println(" ")
    println("decomposition algorithm:")
    println("*" ^ 50)
    println(" algorithm          : " * "Adapt_Bend")
    println(" case               : $(case)" )
    println(" investment  nodes  : $(maximum(Ms.I0))"   )
    println(" operational nodes  : $(maximum(Ms.I ))"    )
    println("-" ^ 50)

end

function string_fcn_ABC_AB(j::Int64,
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

function print_info_AB(J::J_type_AB,
                       B::B_type_AB)

    j = length(B.L)
    Δ = J.Δ
    t = B.t[j]
    str = string_fcn_ABC_AB(j,Δ,t)
    println(str)

end

function str_fcn1_AB(x_ay::Vector{Int64},
                     tech::String)::String

    I = length(x_ay)
    str = " " * tech * " " ^ (9-length(tech))
    for i in 1:I
        x = string(x_ay[i])
        str *= (" " ^ (7-length(x)) * x)
    end

    return str
end

function print_end_summary_AB(Ms::Ms_type,
                               U::U_type,
                               R::R_type_AB,
                               B::B_type_AB,
                            case::Int64)

    println(""); println("*" ^ 75); println("*" ^ 75)
    x0 = zeros(U.nx0,maximum(Ms.I0)); x0 .= getvalue.(R.m[:x0][:,:])
    x1_AB,x2_AB = convert(Array{Int64,1},round.(x0[:,1];digits=0)),convert(Array{Int64,2},round.(x0[:,2:end];digits=0))
    U_AB,L_AB,t_AB,tA_AB,tB_AB,tC_AB = B.U,B.L,B.t,B.tA,B.tC,B.tD
    str_tech = ["coal","coalccs","OCGT","CCGT","diesel","nuclear","pumpL","pumpH","lithium","onwind","offwind","solar"]
    int_A = " tech.    "
    for i in 1:length(x1_AB[1,:])
        st  = "ω$i"
        int_A *= " " ^ (7 - length(st)) * st
    end
    if (case < 4)
        int_B = " tech.    "
        for i in 1:length(x2_AB[1,:]) st  = "ω$i"; int_B *= " " ^ (7 - length(st)) * st; end
        Ni,No = 3^(case -1)+1,3^(2case -2)+3^(case -1)
        obj = round(U_AB[end] / 1e11;digits=3)
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
        for p in 1:12 println(str_fcn1_AB(x1_AB[p,:],str_tech[p])) end
        println(" " * "-" ^ (length(int_A)-1)); println(" ")
        println(" optimal investments @ 5 years"); println(" " * "-" ^ (length(int_B)-1)); println(int_B); println(" " * "-" ^ (length(int_B)-1))
        for p in 1:12 println(str_fcn1_AB(x2_AB[p,:],str_tech[p])) end
        println(" " * "-" ^ (length(int_B)-1)); println(" ")
    end
    if case==4
        int_B = [" tech.    "," tech.    "," tech.    "]
        for j in 1:3 for i in 1:Int64(length(x2_AB[1,:])/3) st  = "ω$(9*(j-1)+i)"; int_B[j] *= " " ^ (7 - length(st)) * st; end end
        Ni,No = 3^(case-1)+1,3^(2case-2)+3^(case-1)
        obj = round(U_AB[end] / 1e11; digits=3)
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
        for p in 1:12 println(str_fcn1_AB(x1_AB[p,:],str_tech[p])) end
        println(" " * "-" ^ (length(int_A)-1)); println(" ")
        println(" optimal investments @ 5 years")
        for j in 1:3
            println(" " * "-" ^ (length(int_B[j])-1)); println(int_B[j]); println(" " * "-" ^ (length(int_B[j])-1))
            for p in 1:12 println(str_fcn1_AB(x2_AB[p,9*(j-1)+1:9*j],str_tech[p])) end
            println(" " * "-" ^ (length(int_B[j])-1)); println(" ")
        end
    end
    println(" computational results:")
    println(" " * "-" ^ 67)
    println(" ϵ-target (%)" * " | " * "iters     time (s)" * "   " * "RMP (%)   SP (%)   Oracles (%)")
    println(" " * "-" ^ 67)
    for ϵ in [1.00,0.10,0.01]
       Δ_AB    = (U_AB.-L_AB)./U_AB*100
       N_AB    = findmin(max.(Δ_AB .- ϵ,0))[2]
       tt_AB,ta_AB,tb_AB,tc_AB = sum(t_AB[1:N_AB]),sum(tA_AB[1:N_AB]),sum(tB_AB[1:N_AB]),sum(tC_AB[1:N_AB])
       RMP_AB,SP_AB,O_AB = round(ta_AB/tt_AB*100;digits=2),round(tb_AB/tt_AB*100;digits=2),round(tc_AB/tt_AB*100;digits=2)
       str_ϵ   = " $(round(ϵ;digits=2))"; str_ϵ   *= "0" ^ (5-length(str_ϵ))
       str_N_AB = "$(N_AB)"  ; str_N_AB *= " " ^ (10-length(str_N_AB))
       str_t_AB = "$(Int64(round(tt_AB;digits=0)))" ; str_t_AB *= " " ^ (9-length(str_t_AB))
       str_R_AB     = "$(RMP_AB)"; str_R_AB = " " ^ (3-findfirst(isequal('.'),str_R_AB)) * str_R_AB; str_R_AB *= "0" ^ (5-length(str_R_AB)); str_R_AB *= " " ^ (10-length(str_R_AB))
       str_S_AB     = "$(SP_AB)" ; str_S_AB = " " ^ (4-findfirst(isequal('.'),str_S_AB)) * str_S_AB; str_S_AB *= "0" ^ (6-length(str_S_AB)); str_S_AB *= " " ^ (9-length(str_S_AB))
       str_O_AB     = "$(O_AB)"  ; str_O_AB = " " ^ (3-findfirst(isequal('.'),str_O_AB)) * str_O_AB; str_O_AB *= "0" ^ (5-length(str_O_AB))

       println(str_ϵ * "         | " * str_N_AB * str_t_AB * "  " * str_R_AB * str_S_AB * str_O_AB)
    end
    println(" " * "-" ^ 67)
    println(""); println("*" ^ 75); println("*" ^ 75)
end
