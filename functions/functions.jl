function run_benders(Ms::Ms_type,
                     Mp::Mp_type,
                     Ps::Ps_type,
                     Pp::Pp_type,
                      U::U_type,
                   case::Int64,
              algorithm::Int64)

    if (algorithm==1)
        R,E,J,B = generate_structs_st(Ms,Mp,Ps,Pp,U)
        print_init_st(case,Ms)
        for it in 1:ITmax
            R,E,J,B = do_step_st(Mp,U,R,E,J,B)
            print_info_st(J,B)
            (J.Δ <= ϵ) ? break : nothing
        end
        print_end_summary_st(Ms,U,R,B,case)
    end

    if (algorithm==2)
        R,E,O,J,B,S,T = generate_structs_or(Ms,Mp,Ps,Pp,U)
        print_init_or(case,Ms)
        E,S,O,J,T = do_step0_or(Ms,Mp,U,E,S,O,J,T)
        for it in 1:ITmax
            R,E,O,S,B,J,T = do_step_or(Mp,U,R,E,O,S,B,J,T)
            print_info_or(J,B)
            (J.Δ <= ϵ) ? break : nothing
        end
        print_end_summary_or(Ms,U,R,B,case)
    end
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
    println(" algorithm 1 -> standard")
    println(" algorithm 2 -> with adaptive oracles")
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
