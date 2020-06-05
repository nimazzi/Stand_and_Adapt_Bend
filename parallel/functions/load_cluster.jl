using Distributed

# */------------------------------------------------------------------------------------------/*

loc_machine = Dict("hostip" => "")

machine_1   = Dict("hostip" => "000.000.000.000",
                   "grblcs" => "path_to_gurobi_license",
                   "jlpath" => "path_to_julia_exec")                  

# */------------------------------------------------------------------------------------------/*

cluster = Dict("machine0" => loc_machine, "machine1"=> machine_1)

# */------------------------------------------------------------------------------------------/*

ipkeys = collect(keys(cluster))
wpools = Dict{String,Array{Int64,1}}()
for k in ipkeys
    nw = get_wrkr(k,cluster[k],cs)
    if cluster[k]["hostip"] == ""
        nw > 0 ? wp = addprocs(nw;topology=:master_worker) : wp = zeros(Int64,0)
    else
        nw > 0 ? wp = addprocs([(cluster[k]["hostip"],nw)];tunnel=true,topology=:master_worker,exename=cluster[k]["jlpath"]) : wp = zeros(Int64,0)
    end
    wpools[k] = wp
end
# */------------------------------------------------------------------------------------------/*

@everywhere using Pkg
@everywhere using JuMP
@everywhere using Suppressor
@everywhere using ParallelDataTransfer
@everywhere using Gurobi

# */------------------------------------------------------------------------------------------/*

for k in ipkeys
    if cluster[k]["hostip"] != ""  
        for w in wpools[k]
            @spawnat(w, Core.eval(Main, Expr(:(=), :grb_lic, cluster[k]["grblcs"])))
            @spawnat w ENV["GRB_LICENSE_FILE"] = grb_lic
            @spawnat w @suppress Pkg.build("Gurobi")
        end
    end
end

# */------------------------------------------------------------------------------------------/*

@everywhere gurobi_env = @suppress Gurobi.Env()
@everywhere Pkg.instantiate()

# */------------------------------------------------------------------------------------------/*
