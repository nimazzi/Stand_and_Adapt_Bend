mutable struct J_type_st
        x::Array{Float64,2}
        β::Vector{Float64}
        θ::Vector{Float64}
        λ::Array{Float64,2}
        Δ::Float64
end

mutable struct B_type_st
    L::Vector{Float64}
    U::Vector{Float64}
    t::Vector{Float64}
   ta::Vector{Float64}
   tb::Vector{Float64}
   tc::Vector{Float64}
end

mutable struct R_type_st
     m::JuMP.Model
  vars::Array{JuMP.Variable,2}
  cons::Array{JuMP.ConstraintRef,2}
end

mutable struct E_type_st
     m::JuMP.Model
  vars::Array{JuMP.Variable,1}
end

function generate_structs_st(Ms::Ms_type,
                             Mp::Mp_type,
                             Ps::Ps_type,
                             Pp::Pp_type,
                              U::U_type)::Tuple{R_type_st,E_type_st,J_type_st,B_type_st}

    m = RMP_model(Ms,Mp,U)
    vars = Array{JuMP.Variable,2}(undef,U.nx+1,U.nI)
    for i in 1:U.nI vars[:,i] = vcat([m[:β][i]],m[:x][:,i]) end
    @constraintref cons[1:ITmax,1:U.nI]
    #   R_type_st(m,vars,cons)
    R = R_type_st(m,vars,cons)
    solve(R.m)

    m = SP_model(Ps,Pp,U)
    vars = vcat(m[:c0],m[:ϕ][:])
    #   E_type_st(m,vars)
    E = E_type_st(m,vars)
    solve(E.m)

    #   J_type_st(x                     ,β          ,θ          ,λ               ,Δ )
    J = J_type_st(getvalue(R.m[:x][:,:]),zeros(U.nI),zeros(U.nI),zeros(U.nx,U.nI),0.)

    #   B_type_st(L       ,U       ,t       ,ta      ,tb      ,tc      )
    B = B_type_st(zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0))

    return R,E,J,B
end

function do_step_A_st(R::R_type_st,
                      J::J_type_st,
                      B::B_type_st)::Tuple{R_type_st,J_type_st,B_type_st}

    solve(R.m)
    J.β .= getvalue(R.m[:β][:])
    J.x .= getvalue(R.m[:x][:,:])
    push!(B.L,getobjectivevalue(R.m))

    return R,J,B
end

function do_step_B_st(U::U_type,
                      E::E_type_st,
                      J::J_type_st)::Tuple{E_type_st,J_type_st}

    for i in 1:U.nI
      @objective(E.m, :Min, AffExpr(E.vars,vcat([1.],U.C[i,:]),.0))
      JuMP.setRHS.(E.m[:λ][:],J.x[:,i])
      solve(E.m)
      J.θ[i] = getobjectivevalue(E.m)
      J.λ[:,i] .= getdual.(E.m[:λ])
    end

    return E,J
end

function do_step_C_st(Mp::Mp_type,
                       U::U_type,
                       J::J_type_st,
                       R::R_type_st,
                       B::B_type_st)::Tuple{R_type_st,B_type_st}

    (length(B.U)==0) ? push!(B.U,getvalue(R.m[:f])+Mp.κ*sum(Mp.πi.*J.θ)) : push!(B.U,min(B.U[end],getvalue(R.m[:f])+Mp.κ*sum(Mp.πi.*J.θ)))
    J.Δ = (B.U[end]-B.L[end])/B.U[end]*100.
    for i in 1:U.nI R.cons[length(B.L),i] = @constraint(R.m, AffExpr(R.vars[:,i],vcat([1.],-J.λ[:,i]),.0) >= J.θ[i]-sum(J.λ[:,i].*J.x[:,i])) end

    return R,B
end

function B_time_st(B::B_type_st,
                   ta::Float64,
                   tb::Float64,
                   tc::Float64)::B_type_st

    push!(B.ta,ta)
    push!(B.tb,   tb)
    push!(B.tc,      tc)
    push!(B.t ,ta+tb+tc)

    return B
end

function do_step_st(Mp::Mp_type,
                     U::U_type,
                     R::R_type_st,
                     E::E_type_st,
                     J::J_type_st,
                     B::B_type_st)::Tuple{R_type_st,E_type_st,J_type_st,B_type_st}

    ta = @elapsed R,J,B = do_step_A_st(R,J,B)
    tb = @elapsed E,J   = do_step_B_st(U,E,J)
    tc = @elapsed R,B   = do_step_C_st(Mp,U,J,R,B)

    B = B_time_st(B,ta,tb,tc)

    return R,E,J,B
end

function print_init_st(case::Int64,
                         Ms::Ms_type)

    println(" ")
    println("decomposition algorithm:")
    println("*" ^ 50)
    println(" algorithm          : " * "Benders standard")
    println(" case               : $(case)" )
    println(" investment  nodes  : $(maximum(Ms.I0))"   )
    println(" operational nodes  : $(maximum(Ms.I ))"    )
    println("-" ^ 50)

end

function string_fcn_ABC_st(j::Int64,
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

function print_info_st(J::J_type_st,
                       B::B_type_st)

    j = length(B.L)
    Δ = J.Δ
    t = B.t[j]
    str = string_fcn_ABC_st(j,Δ,t)
    println(str)

end

function str_fcn1_st(x_ay::Vector{Int64},
                     tech::String)::String

    I = length(x_ay)
    str = " " * tech * " " ^ (9-length(tech))
    for i in 1:I
        x = string(x_ay[i])
        str *= (" " ^ (7-length(x)) * x)
    end

    return str
end

function print_end_summary_st(Ms::Ms_type,
                               U::U_type,
                               R::R_type_st,
                               B::B_type_st,
                            case::Int64)

    println(""); println("*" ^ 75); println("*" ^ 75)
    x0 = zeros(U.nx0,maximum(Ms.I0)); x0 .= getvalue.(R.m[:x0][:,:])
    x1_st,x2_st = convert(Array{Int64,1},round.(x0[:,1];digits=0)),convert(Array{Int64,2},round.(x0[:,2:end];digits=0))
    U_st,L_st,t_st,tA_st,tB_st = B.U,B.L,B.t,B.ta,B.tb
    str_tech = ["coal","coalccs","OCGT","CCGT","diesel","nuclear","pumpL","pumpH","lithium","onwind","offwind","solar"]
    int_A = " tech.    "
    for i in 1:length(x1_st[1,:])
        st  = "ω$i"
        int_A *= " " ^ (7 - length(st)) * st
    end
    if (case < 4)
        int_B = " tech.    "
        for i in 1:length(x2_st[1,:]) st  = "ω$i"; int_B *= " " ^ (7 - length(st)) * st; end
        Ni,No = 3^(case -1)+1,3^(2case -2)+3^(case -1)
        obj = round(U_st[end] / 1e11;digits=3)
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
        for p in 1:12 println(str_fcn1_st(x1_st[p,:],str_tech[p])) end
        println(" " * "-" ^ (length(int_A)-1)); println(" ")
        println(" optimal investments @ 5 years"); println(" " * "-" ^ (length(int_B)-1)); println(int_B); println(" " * "-" ^ (length(int_B)-1))
        for p in 1:12 println(str_fcn1_st(x2_st[p,:],str_tech[p])) end
        println(" " * "-" ^ (length(int_B)-1)); println(" ")
    end
    if case==4
        int_B = [" tech.    "," tech.    "," tech.    "]
        for j in 1:3 for i in 1:Int64(length(x2_st[1,:])/3) st  = "ω$(9*(j-1)+i)"; int_B[j] *= " " ^ (7 - length(st)) * st; end end
        Ni,No = 3^(case-1)+1,3^(2case-2)+3^(case-1)
        obj = round(U_st[end] / 1e11; digits=3)
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
        for p in 1:12 println(str_fcn1_st(x1_st[p,:],str_tech[p])) end
        println(" " * "-" ^ (length(int_A)-1)); println(" ")
        println(" optimal investments @ 5 years")
        for j in 1:3
            println(" " * "-" ^ (length(int_B[j])-1)); println(int_B[j]); println(" " * "-" ^ (length(int_B[j])-1))
            for p in 1:12 println(str_fcn1_st(x2_st[p,9*(j-1)+1:9*j],str_tech[p])) end
            println(" " * "-" ^ (length(int_B[j])-1)); println(" ")
        end
    end
    println(" computational results:")
    println(" " * "-" ^ 53)
    println(" ϵ-target (%)" * " | " * "iters     time (s)" * "   " * "RMP (%)   SP (%)  ")
    println(" " * "-" ^ 53)
    for ϵ in [1.00,0.10,0.01]
        Δ_st = (U_st.-L_st)./U_st*100
        N_st = findmin(max.(Δ_st .- ϵ,0))[2]
        tt_st,ta_st,tb_st = sum(t_st[1:N_st]),sum(tA_st[1:N_st]),sum(tB_st[1:N_st])
        RMPst,SPst     = round(ta_st/tt_st*100;digits=2),round(tb_st/tt_st*100;digits=2)
        str_ϵ   = " $(round(ϵ;digits=2))"                       ; str_ϵ   *= "0" ^ (5-length(str_ϵ))
        str_Nst = "$(N_st)"                                     ; str_Nst *= " " ^ (10-length(str_Nst))
        str_tst = "$(Int64(round(sum(t_st[1:N_st]);digits=0)))" ; str_tst *= " " ^ (8-length(str_tst))
        str_Rst     = "$(RMPst)"; str_Rst = " " ^ (3-findfirst(isequal('.'),str_Rst)) * str_Rst; str_Rst *= "0" ^ (5-length(str_Rst)); str_Rst *= " " ^ (10-length(str_Rst))
        str_Sst     = "$(SPst)" ; str_Sst = " " ^ (4-findfirst(isequal('.'),str_Sst)) * str_Sst; str_Sst *= "0" ^ (6-length(str_Sst)); str_Sst *= " " ^ (6-length(str_Sst))
        println(str_ϵ * "         | " * str_Nst * str_tst * "   " * str_Rst * str_Sst)
    end
    println(" " * "-" ^ 53)
    println(""); println("*" ^ 75); println("*" ^ 75)
end
