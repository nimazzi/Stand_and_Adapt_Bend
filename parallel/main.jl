cd(dirname(@__FILE__))
include("./functions/load_inputfunctions.jl")
J,δ =1000,.01   
cs,al =  get_params()
include("./functions/load_cluster.jl")
include("./functions/load_stuff.jl")

B,S = generate_structures(cs,al,nworkers(),J,δ)

solve_Benders!(B)
println("")
