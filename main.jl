cd(dirname(@__FILE__))                      # set current folder

println("loading packages...")
using JuMP,Gurobi,JLD2,CSV,LinearAlgebra    # load packages

println("loading functions...")
include("functions/data_file.jl")           # load packages data_file.jl
include("functions/functions.jl")           # load packages functions.jl
include("functions/functions_SB.jl")        # load packages functions_SB.jl
include("functions/functions_AB.jl")        # load packages functions_AB.jl
include("functions/opt_models.jl")          # load packages opt_models.jl

const ITmax = 1000                          # set max number of iterations (-)
const ϵ = .01                               # set tolerance gap (%)
case = get_case()::Int64                    # input case study (1, 2, 3 or 4)
algorithm = get_algorithm()::Int64          # input algorithm type (1 or 2)

println("generating datasets...")
Ms,Mp,Ps,Pp,U = generate_data(case)         # generate struct data

println("precompiling code...")
if (algorithm==1)
    Stand_Bend(case,Ms,Mp,Ps,Pp,U)          # run Stand_Bend algorithm
end
if (algorithm==2)
    Adapt_Bend(case,Ms,Mp,Ps,Pp,U)          # run Adapt_Bend algorithm
end

