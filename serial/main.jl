cd(dirname(@__FILE__))
include("./functions/load_inputfunctions.jl")
J,δ =1000,.01   
cs,al,q,w =  get_params()
include("./functions/load_stuff.jl")

B,S = generate_structures(cs,al,q,w,J,δ);
al == 0 ? solve_deterministic!(B,cs) : solve_Benders!(B,S)
println("")



