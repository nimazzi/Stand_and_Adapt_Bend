using Distributed

# */ --------------------------------------------------------------------------------------------- /* #

function get_case()::Int64

    println("")
    println(" case 0 -> 0 uncertain parameters")
    println(" case 1 -> 1 uncertain parameters")
    println(" case 2 -> 2 uncertain parameters")
    println(" case 3 -> 3 uncertain parameters")
    println(" select case study:")
    cs = ""
    cs = readline()
    while (cs ∉ ["0","1","2","3"])
        println(" error: select 0, 1, 2, or 3")
        println(" select case:")
        cs = readline()
    end

    return parse(Int64,cs)
end

function get_algm()::Int64

    println("")
    println(" algorithm 0 -> deterministic equivalent")
    println(" algorithm 1 -> Stand_Bend (Standard Benders)")
    println(" algorithm 2 -> Adapt_Bend (Benders with Adaptive Oracles)")
    println(" algorithm 3 -> Zaker_Bend (Benders, Zakeri et al.)")
    println(" select Benders-type algorithm:")
    al = ""
    al = readline()
    while (al ∉ ["0","1","2","3"])
        println(" error: select 0, 1, 2, or 3")
        println(" select algorithm:")
        al = readline()
    end

    return parse(Int64,al)
end

function get_parw(cs::Int64)::Int64

    cs == 0 ? W = 2   : nothing
    cs == 1 ? W = 12  : nothing
    cs == 2 ? W = 90  : nothing
    cs == 3 ? W = 756 : nothing
    println("")
    println(" number of subproblems to solve each iter j:")
    println(" input w (0 < w < $W)")
    w = ""
    w = readline()
    while (parse(Float64,w) <= 0 || parse(Float64,w) >= W )
        println(" error: select 0 < w < $W")
        println(" input w:")
        w = readline()
    end

    return parse(Float64,w)
end

function get_parσ()::Float64

    println("")
    println(" select parameter σ:")
    println(" input σ (1 < σ <= 100)")
    σ = ""
    σ = readline()
    while (parse(Float64,σ) <= 1 || parse(Float64,σ) > 100 )
        println(" error: select 1 < σ <= 100")
        println(" input σ:")
        σ = readline()
    end

    return parse(Float64,σ)
end

function get_params()::Tuple{Int64,Int64,Float64,Int64}

    c = get_case()
    a = get_algm()
    a == 2 ? w = get_parw(c) : w =  1
    a == 3 ? q = get_parσ( ) : q = .0
    return c,a,q,w
end

# */ --------------------------------------------------------------------------------------------- /* #

