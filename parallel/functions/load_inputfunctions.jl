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
    println(" algorithm 1 -> Stand_Bend (Standard Benders)")
    println(" algorithm 2 -> Adapt_Bend (Benders with Adaptive Oracles)")
    println(" select Benders-type algorithm:")
    al = ""
    al = readline()
    while (al ∉ ["1","2"])
        println(" error: select 1 or 2")
        println(" select algorithm:")
        al = readline()
    end

    return parse(Int64,al)
end

function get_params()::Tuple{Int64,Int64}

    c =  get_case()
    a =  get_algm()

    return c,a
end

# */ --------------------------------------------------------------------------------------------- /* #

function get_wrkr(n::String,d::Dict,case::Int64)::Int64
    
    d["hostip"] == "" ? hip = "local_machine" : hip = d["hostip"]
    case == 0 ? (workers()[1] == 1 ? nw = 2   - 1 : nw = 2   - 1 - nworkers()) : nothing
    case == 1 ? (workers()[1] == 1 ? nw = 12  - 1 : nw = 12  - 1 - nworkers()) : nothing
    case == 2 ? (workers()[1] == 1 ? nw = 90  - 1 : nw = 90  - 1 - nworkers()) : nothing
    case == 3 ? (workers()[1] == 1 ? nw = 756 - 1 : nw = 756 - 1 - nworkers()) : nothing
    println("")
    println(" select number of workers on $(n) ($(hip))")
    println(" select 0 <= w <= $(nw) :")
    w = ""
    w = readline()
    while (parse(Int64,w) < 0 || parse(Int64,w) > nw )
        println(" error ")
        println(" select 0 <= w <= $(nw) :")
        w = readline()
    end
    
    return parse(Int64,w)
end

# */ --------------------------------------------------------------------------------------------- /* #














# */ --------------------------------------------------------------------------------------------- /* #

