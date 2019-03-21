####################################################################################################
###### Algorithm 1: Stand_Bend (Standard Benders ) #################################################
###### iteration (j)
# solve 𝐑𝐌𝐏 and obtain 𝐱ʲ and βʲ                       \* step A  solve master problem            *\
# set Lʲ:= 𝑓(𝐱ʲ)+∑ᵢβʲᵢ                                   \* step B  compute lower bound            *\
# for i ∈ 𝓘                                             \* step C  solve subproblems
#     solve 𝐒𝐏(i) at (xʲᵢ,cᵢ) and obtain θʲᵢ and λʲᵢ                                              *\
# set Uʲ:=𝑓(𝐱ʲ)+∑ᵢθʲᵢ                                    \* step D  compute upper bound            *\
# for i ∈ 𝓘                                             \* step E  update master problem
#     add cut (xʲᵢ,θʲᵢ,λʲᵢ) to 𝐑𝐌𝐏                                                               *\
####################################################################################################

mutable struct J_type_SB
        x::Array{Float64,2}
        β::Vector{Float64}
        θ::Vector{Float64}
        λ::Array{Float64,2}
        Δ::Float64
end

mutable struct B_type_SB
    L::Vector{Float64}
    U::Vector{Float64}
    t::Vector{Float64}
   tA::Vector{Float64}
   tB::Vector{Float64}
   tC::Vector{Float64}
   tD::Vector{Float64}
   tE::Vector{Float64}
end

mutable struct R_type_SB
     m::JuMP.Model
  vars::Array{JuMP.Variable,2}
  cons::Array{JuMP.ConstraintRef,2}
end

mutable struct E_type_SB
     m::JuMP.Model
  vars::Array{JuMP.Variable,1}
end

function gen_structs_SB(Ms::Ms_type,
                        Mp::Mp_type,
                        Ps::Ps_type,
                        Pp::Pp_type,
                         U::U_type)::Tuple{R_type_SB,E_type_SB,J_type_SB,B_type_SB}

    m = RMP_model(Ms,Mp,U)
    vars = Array{JuMP.Variable,2}(undef,U.nx+1,U.nI)
    for i in 1:U.nI vars[:,i] = vcat([m[:β][i]],m[:x][:,i]) end
    @constraintref cons[1:ITmax,1:U.nI]
    #   R_type_SB(m,vars,cons)
    R = R_type_SB(m,vars,cons)
    solve(R.m)

    m = SP_model(Ps,Pp,U)
    vars = vcat(m[:c0],m[:ϕ][:])
    #   E_type_SB(m,vars)
    E = E_type_SB(m,vars)
    solve(E.m)

    #   J_type_SB(x                     ,β          ,θ          ,λ               ,Δ )
    J = J_type_SB(getvalue(R.m[:x][:,:]),zeros(U.nI),zeros(U.nI),zeros(U.nx,U.nI),0.)

    #   B_type_SB(L       ,U       ,t       ,tA      ,tB      ,tC      ,tD      ,tE      )
    B = B_type_SB(zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0))

    return R,E,J,B
end

function Stand_Bend_step_A(R::R_type_SB,
                           J::J_type_SB)::Tuple{R_type_SB,J_type_SB}

    solve(R.m)
    J.β .= getvalue(R.m[:β][:])
    J.x .= getvalue(R.m[:x][:,:])

    return R,J
end

function Stand_Bend_step_B(R::R_type_SB,
                           B::B_type_SB)::B_type_SB

    push!(B.L,getobjectivevalue(R.m))

    return B
end

function Stand_Bend_step_C(U::U_type,
                           E::E_type_SB,
                           J::J_type_SB)::Tuple{E_type_SB,J_type_SB}

    for i in 1:U.nI
      @objective(E.m, :Min, AffExpr(E.vars,vcat([1.],U.C[i,:]),.0))
      JuMP.setRHS.(E.m[:λ][:],J.x[:,i])
      solve(E.m)
      J.θ[i] = getobjectivevalue(E.m)
      J.λ[:,i] .= getdual.(E.m[:λ])
    end

    return E,J
end

function Stand_Bend_step_D(Mp::Mp_type,
                            U::U_type,
                            J::J_type_SB,
                            R::R_type_SB,
                            B::B_type_SB)::B_type_SB

    (length(B.U)==0) ? push!(B.U,getvalue(R.m[:f])+Mp.κ*sum(Mp.πi.*J.θ)) : push!(B.U,min(B.U[end],getvalue(R.m[:f])+Mp.κ*sum(Mp.πi.*J.θ)))

    return B
end

function Stand_Bend_step_E(U::U_type,
                           J::J_type_SB,
                           R::R_type_SB,
                           B::B_type_SB)::Tuple{R_type_SB,J_type_SB}

    for i in 1:U.nI R.cons[length(B.L),i] = @constraint(R.m, AffExpr(R.vars[:,i],vcat([1.],-J.λ[:,i]),.0) >= J.θ[i]-sum(J.λ[:,i].*J.x[:,i])) end
    J.Δ = (B.U[end]-B.L[end])/B.U[end]*100.

    return R,J
end

function B_time_SB(B::B_type_SB,
                  tA::Float64,
                  tB::Float64,
                  tC::Float64,
                  tD::Float64,
                  tE::Float64)::B_type_SB

    push!(B.tA,tA)
    push!(B.tB,   tB)
    push!(B.tC,      tC)
    push!(B.tD,         tD)
    push!(B.tE,            tE)
    push!(B.t ,tA+tB+tC+tD+tE)

    return B
end

function iter_SB(Mp::Mp_type,
                  U::U_type,
                  R::R_type_SB,
                  E::E_type_SB,
                  J::J_type_SB,
                  B::B_type_SB)::Tuple{R_type_SB,E_type_SB,J_type_SB,B_type_SB}

    tA = @elapsed R,J   = Stand_Bend_step_A(R,J)
    tB = @elapsed B     = Stand_Bend_step_B(R,B)
    tC = @elapsed E,J   = Stand_Bend_step_C(U,E,J)
    tD = @elapsed B     = Stand_Bend_step_D(Mp,U,J,R,B)
    tE = @elapsed R,J   = Stand_Bend_step_E(U,J,R,B)

    B  = B_time_SB(B,tA,tB,tC,tD,tE)

    return R,E,J,B
end

function print_init_SB(case::Int64,
                         Ms::Ms_type)

    println(" ")
    println("decomposition algorithm:")
    println("*" ^ 50)
    println(" algorithm          : " * "Stand_Bend")
    println(" case               : $(case)" )
    println(" investment  nodes  : $(maximum(Ms.I0))"   )
    println(" operational nodes  : $(maximum(Ms.I ))"    )
    println("-" ^ 50)

end

function string_fcn_ABC_SB(j::Int64,
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

function print_info_SB(J::J_type_SB,
                       B::B_type_SB)

    j = length(B.L)
    Δ = J.Δ
    t = B.t[j]
    str = string_fcn_ABC_SB(j,Δ,t)
    println(str)

end

function str_fcn1_SB(x_ay::Vector{Int64},
                     tech::String)::String

    I = length(x_ay)
    str = " " * tech * " " ^ (9-length(tech))
    for i in 1:I
        x = string(x_ay[i])
        str *= (" " ^ (7-length(x)) * x)
    end

    return str
end

function print_end_summary_SB(Ms::Ms_type,
                               U::U_type,
                               R::R_type_SB,
                               B::B_type_SB,
                            case::Int64)

    println(""); println("*" ^ 75); println("*" ^ 75)
    x0 = zeros(U.nx0,maximum(Ms.I0)); x0 .= getvalue.(R.m[:x0][:,:])
    x1_SB,x2_SB = convert(Array{Int64,1},round.(x0[:,1];digits=0)),convert(Array{Int64,2},round.(x0[:,2:end];digits=0))
    U_SB,L_SB,t_SB,tA_SB,tC_SB = B.U,B.L,B.t,B.tA,B.tC
    str_tech = ["coal","coalccs","OCGT","CCGT","diesel","nuclear","pumpL","pumpH","lithium","onwind","offwind","solar"]
    int_A = " tech.    "
    for i in 1:length(x1_SB[1,:])
        st  = "ω$i"
        int_A *= " " ^ (7 - length(st)) * st
    end
    if (case < 4)
        int_B = " tech.    "
        for i in 1:length(x2_SB[1,:]) st  = "ω$i"; int_B *= " " ^ (7 - length(st)) * st; end
        Ni,No = 3^(case -1)+1,3^(2case -2)+3^(case -1)
        obj = round(U_SB[end] / 1e11;digits=3)
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
        for p in 1:12 println(str_fcn1_SB(x1_SB[p,:],str_tech[p])) end
        println(" " * "-" ^ (length(int_A)-1)); println(" ")
        println(" optimal investments @ 5 years"); println(" " * "-" ^ (length(int_B)-1)); println(int_B); println(" " * "-" ^ (length(int_B)-1))
        for p in 1:12 println(str_fcn1_SB(x2_SB[p,:],str_tech[p])) end
        println(" " * "-" ^ (length(int_B)-1)); println(" ")
    end
    if case==4
        int_B = [" tech.    "," tech.    "," tech.    "]
        for j in 1:3 for i in 1:Int64(length(x2_SB[1,:])/3) st  = "ω$(9*(j-1)+i)"; int_B[j] *= " " ^ (7 - length(st)) * st; end end
        Ni,No = 3^(case-1)+1,3^(2case-2)+3^(case-1)
        obj = round(U_SB[end] / 1e11; digits=3)
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
        for p in 1:12 println(str_fcn1_SB(x1_SB[p,:],str_tech[p])) end
        println(" " * "-" ^ (length(int_A)-1)); println(" ")
        println(" optimal investments @ 5 years")
        for j in 1:3
            println(" " * "-" ^ (length(int_B[j])-1)); println(int_B[j]); println(" " * "-" ^ (length(int_B[j])-1))
            for p in 1:12 println(str_fcn1_SB(x2_SB[p,9*(j-1)+1:9*j],str_tech[p])) end
            println(" " * "-" ^ (length(int_B[j])-1)); println(" ")
        end
    end
    println(" computational results:")
    println(" " * "-" ^ 53)
    println(" ϵ-target (%)" * " | " * "iters     time (s)" * "   " * "RMP (%)   SP (%)  ")
    println(" " * "-" ^ 53)
    for ϵ in [1.00,0.10,0.01]
        Δ_SB = (U_SB.-L_SB)./U_SB*100
        N_SB = findmin(max.(Δ_SB .- ϵ,0))[2]
        tt_SB,ta_SB,tc_SB = sum(t_SB[1:N_SB]),sum(tA_SB[1:N_SB]),sum(tC_SB[1:N_SB])
        RMP_SB,SP_SB = round(ta_SB/tt_SB*100;digits=2),round(tc_SB/tt_SB*100;digits=2)
        str_ϵ   = " $(round(ϵ;digits=2))"                       ; str_ϵ   *= "0" ^ (5-length(str_ϵ))
        str_N_SB = "$(N_SB)"                                     ; str_N_SB *= " " ^ (10-length(str_N_SB))
        str_t_SB = "$(Int64(round(sum(t_SB[1:N_SB]);digits=0)))" ; str_t_SB *= " " ^ (8-length(str_t_SB))
        str_R_SB     = "$(RMP_SB)"; str_R_SB = " " ^ (3-findfirst(isequal('.'),str_R_SB)) * str_R_SB; str_R_SB *= "0" ^ (5-length(str_R_SB)); str_R_SB *= " " ^ (10-length(str_R_SB))
        str_S_SB     = "$(SP_SB)" ; str_S_SB = " " ^ (4-findfirst(isequal('.'),str_S_SB)) * str_S_SB; str_S_SB *= "0" ^ (6-length(str_S_SB)); str_S_SB *= " " ^ (6-length(str_S_SB))
        println(str_ϵ * "         | " * str_N_SB * str_t_SB * "   " * str_R_SB * str_S_SB)
    end
    println(" " * "-" ^ 53)
    println(""); println("*" ^ 75); println("*" ^ 75)
end
