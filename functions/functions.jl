function Stand_Bend(case::Int64,
                      Ms::Ms_type,
                      Mp::Mp_type,
                      Ps::Ps_type,
                      Pp::Pp_type,
                       U::U_type)

    R,E,J,B = gen_structs_SB(Ms,Mp,Ps,Pp,U)
    print_init_SB(case,Ms)
    for it in 1:ITmax
        R,E,J,B = iter_SB(Mp,U,R,E,J,B)
        print_info_SB(J,B)
        (J.Δ <= ϵ) ? break : nothing
    end
    print_end_summary_SB(Ms,U,R,B,case)

end

function Adapt_Bend(case::Int64,
                      Ms::Ms_type,
                      Mp::Mp_type,
                      Ps::Ps_type,
                      Pp::Pp_type,
                       U::U_type)

    R,E,O,J,B,S,T = gen_structs_AB(Ms,Mp,Ps,Pp,U)
    print_init_AB(case,Ms)
    E,S,O,J,T = Adapt_Bend_step_0(Ms,Mp,U,E,S,O,J,T)
    for it in 1:ITmax
       R,E,O,S,B,J,T = iter_AB(Mp,U,R,E,O,S,B,J,T)
       print_info_AB(J,B)
       (J.Δ <= ϵ) ? break : nothing
    end
    print_end_summary_AB(Ms,U,R,B,case)

end

function get_case()::Int64

    println("")
    println(" case 1 -> 0 uncertain parameters")
    println(" case 2 -> 1 uncertain parameters")
    println(" case 3 -> 2 uncertain parameters")
    println(" case 4 -> 3 uncertain parameters")
    println(" select case study:")
    cs = ""
    cs = readline()
    while (cs != "1" && cs != "2" && cs != "3" && cs != "4")
        println(" error: select 1, 2, 3, or 4")
        println(" select case:")
        cs = readline()
    end

    return parse(Int64,cs)
end

function get_algorithm()::Int64

    println("")
    println(" algorithm 1 -> Stand_Bend (Standard Benders)")
    println(" algorithm 2 -> Adapt_Bend (Benders with Adaptive Oracles)")
    println(" select Benders-type algorithm:")
    al = ""
    al = readline()
    while (al != "1" && al != "2")
        println(" error: select 1 or 2")
        println(" select algorithm:")
        al = readline()
    end
    println("")

    return parse(Int64,al)
end
