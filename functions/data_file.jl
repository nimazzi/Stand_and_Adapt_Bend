# define structures
struct Ms_type # master problem sets
        P::UnitRange{Int64}       # set of technologies
       I0::UnitRange{Int64}       # set of "investment" nodes
        I::UnitRange{Int64}       # set of "operational" nodes
end

struct Mp_type # master problem parameters
        κ::Float64                # years of operational problem
   x_hist::Array{Float64,2}       # historical installed capacity (GW)
    x_max::Vector{Float64}        # maximum accumulated capacity (GW)
    c_inv::Array{Float64,2}       # investment cost (£/GW)
    c_fix::Vector{Float64}        # fix OM cost (£/GWyr)
      πi0::Vector{Float64}        # probability associated to "investment" node i0 (-)
       πi::Vector{Float64}        # probability associated to "operational" node i (-)
    itoi0::Vector{Array{Int64,1}} # map "operational" node i → "investment" node i0
end

struct Ps_type   # subproblem sets
        P::UnitRange{Int64}       # set of technologies
        G::UnitRange{Int64}       # set of conventional technologies
        B::UnitRange{Int64}       # set of storage technologies
        R::UnitRange{Int64}       # set of renewable technologies
        S::UnitRange{Int64}       # set of seasons
        H::UnitRange{Int64}       # hours within a season
       b0::Int64
       r0::Int64
end

struct Pp_type # subproblem parameters
 c_OMvarG::Vector{Float64}        # var OM cost (conv. tech.) (£/MWh)
  c_fuelG::Vector{Float64}        # fuel cost (conv. tech.) (£/MWh)
  em_co2G::Vector{Float64}        # CO2 emissions (conv. tech.) (tCO₂/MWh)
       ηG::Vector{Float64}        # efficiency (conv. tech.) (-)
    rampG::Vector{Float64}        # ramping limitations (conv. tech.) (MW)
       ηB::Vector{Float64}        # efficiency (storage tech.) (-)
       PB::Vector{Float64}        # ch/disch power (storage tech.) (MW/MW)
       PR::Array{Float64,3}       # renewable energy prod (re. tech.) (MW/MW)
       PD::Array{Float64,2}       # energy demand (MW)
       νD::Float64                # scaling demand (-)
   c_shed::Float64                # shedding demand cost (£/MW)
        α::Float64                # seasonal weight (-)
    c_co2::Float64                # CO2 cost (conv. tech.) (£/tCO₂)
  co2_lim::Float64                # co2 yearly limit (tCO₂)
     ρco2::Float64                # scaling co2 limit (-)
end

struct U_type
     nI::Int64
     nx::Int64
    nx0::Int64
     nh::Int64
     nc::Int64
      H::Array{Float64,2}
      C::Array{Float64,2}
end

function generate_data(case::Int64)::Tuple{Ms_type,Mp_type,Ps_type,Pp_type,U_type}

    df_sets = CSV.read("$(pwd())/input/df_sets.csv" ;delim=",",use_mmap=false,nullable=false)
    P  =               1:df_sets[:P ][1]
    G  = df_sets[:G0][1]:df_sets[:G1][1]
    B  = df_sets[:B0][1]:df_sets[:B1][1]
    R  = df_sets[:R0][1]:df_sets[:R1][1]
    S  =               1:df_sets[:S ][1]
    Hs =               1:df_sets[:Hs][1]
    b0 =                 df_sets[:B0][1]-1
    r0 =                 df_sets[:R0][1]-1

    df_unc_sets = CSV.read("$(pwd())/input/df_unc_sets_case$(case).csv" ;delim=",",use_mmap=false,nullable=false)
    NI1 = df_unc_sets[:NΩ1][1]
    NI2 = df_unc_sets[:NΩ2][1]
    NI3 = df_unc_sets[:NΩ3][1]
    NI0 = NI1 + NI2
    NI  = NI2 + NI3

    df_unc_params = CSV.read("$(pwd())/input/df_unc_params_case$(case).csv" ;delim=",",use_mmap=false,nullable=false)
    H      = zeros(NI,2)
    C      = zeros(NI,2)
    H[:,1] = df_unc_params[:ρco2 ]
    H[:,2] = df_unc_params[:νD   ]
    C[:,1] = df_unc_params[:c_co2]
    C[:,2] = df_unc_params[:c_ur ]

    df_inv_params_P = CSV.read("$(pwd())/input/df_inv_params_P.csv" ;delim=",",use_mmap=false,nullable=false)
    c_inv_0   = df_inv_params_P[:c_inv_0  ][:]
    c_inv_5   = df_inv_params_P[:c_inv_5  ][:]
    c_fix     = df_inv_params_P[:c_fix    ][:].*10^3
    x_hist_5  = df_inv_params_P[:x_hist_5 ][:]
    x_hist_10 = df_inv_params_P[:x_hist_10][:]
    x_max     = df_inv_params_P[:x_max    ][:].*10^-3
    c_inv     = hcat(repeat(c_inv_0 ,1,NI1),repeat(c_inv_5  ,1,NI2)).*10^3
    x_hist    = hcat(repeat(x_hist_5,1,NI2),repeat(x_hist_10,1,NI3)).*10^-3
    πi0       = vcat(1/NI1*ones(NI1),1/NI2*ones(NI2))
    πi        = vcat(1/NI2*ones(NI2),1/NI3*ones(NI3))
    itoi0     = Vector{Array{Int64,1}}()  # [ωM,..] = ΩStoM[ωS]
    for i in 1:NI2 push!(itoi0,[1]) end
    for i in (NI2+1):(NI2+NI3) push!(itoi0,[1,Int64(((i-NI2-1)-(i-NI2-1)%NI2)/NI2+1)+1]) end

    df_inv_params_O = CSV.read("$(pwd())/input/df_inv_params_O.csv" ;delim=",",use_mmap=false,nullable=false)
    κ = df_inv_params_O[:κ][1]

    df_oper_params_G = CSV.read("$(pwd())/input/df_oper_params_G.csv" ;delim=",",use_mmap=false,nullable=false)
    c_OMvarG = df_oper_params_G[:c_varOM][:]
    c_fuelG  = df_oper_params_G[:c_fuel ][:]
    em_co2G  = df_oper_params_G[:em_co2 ][:]
    ηG       = df_oper_params_G[:η      ][:]
    rampG    = df_oper_params_G[:ramp   ][:]

    df_oper_params_B = CSV.read("$(pwd())/input/df_oper_params_B.csv" ;delim=",",use_mmap=false,nullable=false)
    ηB = df_oper_params_B[:η][:]
    PB = df_oper_params_B[:P][:]

    df_oper_params_R = CSV.read("$(pwd())/input/df_oper_params_R.csv" ;delim=",",use_mmap=false,nullable=false)
    PR = zeros(maximum(R)-r0+1,maximum(S),maximum(Hs))
    PR[1,1,:] = df_oper_params_R[:P1_S1]
    PR[1,2,:] = df_oper_params_R[:P1_S2]
    PR[1,3,:] = df_oper_params_R[:P1_S3]
    PR[1,4,:] = df_oper_params_R[:P1_S4]
    PR[2,1,:] = df_oper_params_R[:P2_S1]
    PR[2,2,:] = df_oper_params_R[:P2_S2]
    PR[2,3,:] = df_oper_params_R[:P2_S3]
    PR[2,4,:] = df_oper_params_R[:P2_S4]
    PR[3,1,:] = df_oper_params_R[:P3_S1]
    PR[3,2,:] = df_oper_params_R[:P3_S2]
    PR[3,3,:] = df_oper_params_R[:P3_S3]
    PR[3,4,:] = df_oper_params_R[:P3_S4]

    df_oper_params_D = CSV.read("$(pwd())/input/df_oper_params_D.csv" ;delim=",",use_mmap=false,nullable=false)
    PD = zeros(maximum(S),maximum(Hs))
    PD[1,:] = df_oper_params_D[:P_S1]
    PD[2,:] = df_oper_params_D[:P_S2]
    PD[3,:] = df_oper_params_D[:P_S3]
    PD[4,:] = df_oper_params_D[:P_S4]

    df_oper_params_O = CSV.read("$(pwd())/input/df_oper_params_O.csv" ;delim=",",use_mmap=false,nullable=false)
    νD      = df_oper_params_O[:νD     ][1]
    c_shed  = df_oper_params_O[:c_shed ][1]
    α       = df_oper_params_O[:α      ][1]
    c_co2   = df_oper_params_O[:c_co2  ][1]
    co2_lim = df_oper_params_O[:co2_lim][1]
    ρco2    = df_oper_params_O[:ρco2   ][1]

    #   U_type(nI,nx          ,nx0       ,nh,nc,H,C)
    U = U_type(NI,maximum(P)+2,maximum(P),2 , 2,H,C)

    #    Ms_type(P,I0   ,I   )
    Ms = Ms_type(P,1:NI0,1:NI)

    #    Mp_type(κ,x_hist,x_max,c_inv,c_fix,πi0,πi,itoi0)
    Mp = Mp_type(κ,x_hist,x_max,c_inv,c_fix,πi0,πi,itoi0)

    #    Ps_type(P,G,B,R,S,Hs,b0,r0)
    Ps = Ps_type(P,G,B,R,S,Hs,b0,r0)

    #    Pp_type(c_OMvarG,c_fuelG,em_co2G,ηG,rampG,ηB,PB,PR,PD,νD,c_shed,α,c_co2,co2_lim,ρco2)
    Pp = Pp_type(c_OMvarG,c_fuelG,em_co2G,ηG,rampG,ηB,PB,PR,PD,νD,c_shed,α,c_co2,co2_lim,ρco2)

    return Ms,Mp,Ps,Pp,U
end
