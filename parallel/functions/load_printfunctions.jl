# */ --------------------------------------------------------------------------------------------- /* #

function print0(b::B1_type)
    
    println(" ")
    println(" */"*"-"^ 32*"/*")
    println(" algorithm          : " * "stand_bend")
    println(" case               : $(b.data.case)")
    println(" investment  nodes  : $(b.data.ms.I0[end])")
    println(" operational nodes  : $(b.data.ms.I[end])")
    println(" workers            : $(b.data.wrks)")
    println(" */"*"-"^ 32*"/*")
    
end

function print0(b::B2_type)
    
    println(" ")
    println(" */"*"-"^ 32*"/*")
    println(" algorithm          : " * "adapt_bend")
    println(" case               : $(b.data.case)")
    println(" investment  nodes  : $(b.data.ms.I0[end])")
    println(" operational nodes  : $(b.data.ms.I[end])")
    println(" workers            : $(b.data.wrks)")
    println(" */"*"-"^ 32*"/*")
    
end

# */ --------------------------------------------------------------------------------------------- /* #

function print_info(b::B1_type)

    a1 = "$(round(b.temp.Δ;digits=3))" * "0" ^ (5 - length("$(round(b.temp.Δ%1;digits=3))"))
    a2 = "$(round(sum(b.hist.T[b.hist.k,:]);digits=2))"
    str = " k =" * (" " ^ (4-length("$(b.hist.k)"))) * "$(b.hist.k), "
    str *= "δ =" * (" " ^ (7-length(a1))) * "$(a1) %, "
    str *= "t =" * (" " ^ (8-length(a2))) * "$(a2) s"
    println(str)

end

function print_info(b::B2_type)

    a1 = "$(round(b.temp.Δ;digits=3))" * "0" ^ (5 - length("$(round(b.temp.Δ%1;digits=3))"))
    a2 = "$(round(sum(b.hist.T[b.hist.k,:]);digits=2))"
    str = " k =" * (" " ^ (4-length("$(b.hist.k)"))) * "$(b.hist.k), "
    str *= "δ =" * (" " ^ (7-length(a1))) * "$(a1) %, "
    str *= "t =" * (" " ^ (8-length(a2))) * "$(a2) s"
    println(str)

end

# */ --------------------------------------------------------------------------------------------- /* #

function str_fcn1(x::Array{Float64,1},s::String)::String

    str = " "*s*" "^(9-length(s))
    for i in 1:length(x)
        str *= (" "^(7-length(string(x[i])))*string(x[i]))
    end

    return str
end

function print_summary(b)

    print_summary_a(b)
    print_summary_b(b)

end

function print_summary_a(b)

    println(" */"*"-"^ 32*"/*")
    println(" ")
    println(" */"*"-"^ 68*"/*")
    x0 = value.(b.rmp.m[:x0])
    x1 = round.(x0[:,1];digits=1)
    x2 = round.(x0[:,2:end];digits=1)
    str_tech = ["coal","coalccs","OCGT","CCGT","diesel","nuclear","pumpL","pumpH","lithium","onwind","offwind","solar"]
    int1 = " tech."*" "^9*"i1"
    println(" ")
    println(" co2 emission limit : $( b.data.case <= 1 ? "known" : "uncertain" )")
    println(" co2 emission cost  : $( b.data.case <= 2 ? "known" : "uncertain" )")
    println(" uranium cost       : $( b.data.case <= 3 ? "known" : "uncertain" )")
    println(" investment nodes   : $( b.data.ms.I0[end] )")
    println(" operational nodes  : $( b.data.ms.I[ end] )")
    println(" optimal objective  : $( round(b.hist.U[b.hist.k]*exp10(-11);digits=3) ) x 10^11 £")
    println(" ")
    println(" */"*"-"^ 68*"/*")
    println(" ")
    println(" optimal investments (GW) @ 0 years")
    println(" " * "-" ^ (length(int1)-1))
    println(int1)
    println(" " * "-" ^ (length(int1)-1))
    for p in b.data.ms.P 
        println(str_fcn1(x1[p,:],str_tech[p])) 
    end
    println(" "*"-"^(length(int1)-1))
    println(" ")
    println(" optimal investments (GW) @ 5 years")
    if (b.data.case<=3)
        int2 = " tech.    "
        for n in 1:size(x2)[2] 
            int2 *= " "^(7-length("i$n"))*"i$n"
        end
        println(" " * "-" ^ (length(int2)-1))
        println(int2)
        println(" " * "-" ^ (length(int2)-1))
        for p in b.data.ms.P 
            println(str_fcn1(x2[p,:],str_tech[p])) 
        end
        println(" "*"-"^(length(int2)-1))
    else
        int2 = [" tech.    "," tech.    "," tech.    "]
        for j in 1:3 
            for i in 1:Int64(size(x2)[2]/3) 
                int2[j] *= " "^(7-length("i$(9*(j-1)+i)"))*"i$(9*(j-1)+i)"
            end 
        end
        for j in 1:3
            println(" " * "-" ^ (length(int2[j])-1))
            println(int2[j])
            println(" " * "-" ^ (length(int2[j])-1))
            for p in b.data.ms.P 
                println(str_fcn1(x2[p,9*(j-1)+1:9*j],str_tech[p])) 
            end
            println(" " * "-" ^ (length(int2[j])-1))
        end
    end
    println(" ")
    println(" */"*"-"^ 68*"/*")
end

function print_summary_b(b)

    println(" ")
    println(" computational results:")
    println(" " * "-" ^ 67)
    println(" ϵ-target (%)" * " | " * "iters     time (s)" * "   " * "RMP (%)   SP (%)   Oracles (%)")
    println(" " * "-" ^ 67)
    δ = exp10(2)*(b.hist.U[1:b.hist.k].-b.hist.L[1:b.hist.k])./b.hist.U[1:b.hist.k]
    for e in [1.00, 0.10, 0.01]
        n =findmin(max.(δ.-e,0))[2]
        t  = sum(b.hist.T[1:n,:])
        tm = sum(b.hist.T[1:n,1])
        ts = sum(b.hist.T[1:n,3])
        to = sum(b.hist.T[1:n,4])
        rm = round(exp10(2)*tm/t;digits=2)
        rs = round(exp10(2)*ts/t;digits=2)
        ro = round(exp10(2)*to/t;digits=2)
        str_e = " $(round(e;digits=2))"
        str_n = "$(n)"
        str_t = "$(Int64(round(t;digits=0)))"
        str_e *= "0" ^ (5-length(str_e))
        str_n *= " " ^ (10-length(str_n))
        str_t *= " " ^ (9-length(str_t))
        str_m = "$(rm)"
        str_s = "$(rs)"
        str_o = "$(ro)"
        str_m = " " ^ (3-findfirst(isequal('.'),str_m)) * str_m
        str_s = " " ^ (4-findfirst(isequal('.'),str_s)) * str_s
        str_o = " " ^ (3-findfirst(isequal('.'),str_o)) * str_o
        str_m *= "0" ^ (5-length(str_m))
        str_s *= "0" ^ (6-length(str_s))
        str_o *= "0" ^ (5-length(str_o))
        str_m *= " " ^ (10-length(str_m))
        str_s *= " " ^ (9-length(str_s))
        println(str_e * "         | " * str_n * str_t * "  " * str_m * str_s * str_o)
    end
    println(" " * "-" ^ 67)
    println(" ")
    println(" */"*"-"^ 68*"/*")

end
