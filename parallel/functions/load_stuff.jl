@everywhere push!(LOAD_PATH,"$(pwd())/structures")
@everywhere push!(LOAD_PATH,"$(pwd())/functions")
include("./load_packages.jl")
include("./load_structures.jl")
include("./load_data.jl")
include("./load_functions.jl")
include("./load_printfunctions.jl")