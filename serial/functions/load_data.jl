function load_data(case::Int64)::Tuple{ms_type,mp_type,ps_type,pp_type,u_type}

    df_sets = CSV.read("$(pwd())/data/df_sets.csv" ;delim=",",use_mmap=false)
    P  = 1:df_sets[1,:P]                                             # set of technologies 𝒫
    G  = 1:df_sets[1,:G]                                             # set of conventional technologies 𝒢
    B  = 1:df_sets[1,:B]                                             # set of storage technologies ℬ
    R  = 1:df_sets[1,:R]                                             # set of renewable technologies ℛ
    S  = 1:df_sets[1,:S]                                             # set of seasons 𝒮
    Hs = 1:df_sets[1,:H]                                             # set of hours ℋ per seasons
    G0 = 1:df_sets[1,:G0]                                            # set of conventional technologies without nuclear
    gn = df_sets[1,:gn]                                              # index of nuclear generation
    b0 = df_sets[1,:b0]                                              # starting index of storage technologies
    r0 = df_sets[1,:r0]                                              # starting index of renewable technologies

    df_unc_sets = CSV.read("$(pwd())/data/df_unc_sets_case$(case).csv" ;delim=",",use_mmap=false)
    NI1 = df_unc_sets[1,:NΩ1]                                        # number of stochastic nodes at present time
    NI2 = df_unc_sets[1,:NΩ2]                                        # number of stochastic nodes in 5 years
    NI3 = df_unc_sets[1,:NΩ3]                                        # number of stochastic nodes in 10 years
    NI0 = NI1 + NI2                                                  # number of investiment nodes 
    NI  = NI2 + NI3                                                  # number of operational nodes 

    df_unc_params = CSV.read("$(pwd())/data/df_unc_params_case$(case).csv" ;delim=",",use_mmap=false)
    H      = zeros(NI,2)                                             # create empty matrix of rhs terms uncertain paramters H
    C      = zeros(NI,2)                                             # create empty matrix of cost coefficient uncertain paramters C
    H[:,1] = df_unc_params[:,:ρco2 ]                                 # set values of CO₂ budget (H)
    H[:,2] = -df_unc_params[:,:νD  ]                                 # set values of demand growth (H)
    C[:,1] = df_unc_params[:,:c_co2]                                 # set values of CO₂ emission cost (C)
    C[:,2] = df_unc_params[:,:c_ur ]                                 # set values of nuclear fuel cost (C)

    df_inv_params_P = CSV.read("$(pwd())/data/df_inv_params_P.csv" ;delim=",",use_mmap=false)
    c_inv_0   = df_inv_params_P[:,:c_inv_0].*exp10(3)                # investment cost of technologies at present time (£/GW)
    c_inv_5   = df_inv_params_P[:,:c_inv_5  ].*exp10(3)              # investment cost of technologies in 5 years (£/GW)
    c_fix     = df_inv_params_P[:,:c_fix    ].*exp10(3)              # fixed O&M cost of technologies (£/GWyr)
    x_hist_5  = df_inv_params_P[:,:x_hist_5 ].*exp10(-3)             # historical values of technologies in 5 years (GW)
    x_hist_10 = df_inv_params_P[:,:x_hist_10].*exp10(-3)             # historical values of technologies in 10 years (GW)
    x_max     = df_inv_params_P[:,:x_max    ].*exp10(-3)             # maximum capacity of technologies (GW)
    c_inv     = hcat(repeat(c_inv_0 ,1,NI1),repeat(c_inv_5  ,1,NI2)) # investment cost of technologies (£/GW)
    x_hist    = hcat(repeat(x_hist_5,1,NI2),repeat(x_hist_10,1,NI3)) # historical values of technologies (GW)
    πi0       = vcat(1/NI1*ones(NI1),1/NI2*ones(NI2))                # probabilities associated to investiment nodes 
    πi        = vcat(1/NI2*ones(NI2),1/NI3*ones(NI3))                # probabilities associated to operational nodes 
    itoi0     = Vector{Array{Int64,1}}()                             # map that, for given operational node i, it gives its ancestor investiment nodes
    for i in 1:NI2 push!(itoi0,[1]) end
    for i in (NI2+1):(NI2+NI3) push!(itoi0,[1,Int64(((i-NI2-1)-(i-NI2-1)%NI2)/NI2+1)+1]) end

    df_inv_params_O = CSV.read("$(pwd())/data/df_inv_params_O.csv" ;delim=",",use_mmap=false)
    κ = df_inv_params_O[1,:κ]                                        # years of operational problem (yr)

    df_oper_params_G = CSV.read("$(pwd())/data/df_oper_params_G.csv" ;delim=",",use_mmap=false)
    c_OMvarG = df_oper_params_G[:,:c_varOM]                          # set variable OM cost (conventional technologies) (£/MWh)        
    c_fuelG  = df_oper_params_G[:,:c_fuel ]                          # set fuel cost (conventional technologies) (£/MWh)
    em_co2G  = df_oper_params_G[:,:em_co2 ]                          # set CO2 emissions (conventional technologies) (tCO₂/MWh)
    ηG       = df_oper_params_G[:,:η      ]                          # set efficiency (conventional technologies)
    rampG    = df_oper_params_G[:,:ramp   ]                          # ramping limitations (conventional technologies) (MW)

    df_oper_params_B = CSV.read("$(pwd())/data/df_oper_params_B.csv" ;delim=",",use_mmap=false)
    ηB = df_oper_params_B[:,:η]                                      # set efficiency (storage technologies) 
    PB = df_oper_params_B[:,:P]                                      # set ch/disch power (storage technologies) (MW/MW)

    df_oper_params_R = CSV.read("$(pwd())/data/df_oper_params_R.csv" ;delim=",",use_mmap=false)
    PR = zeros(R[end],S[end],Hs[end])                                # create empty array of renewable energy production (MW/MW) 
    PR[1,1,:] = df_oper_params_R[:,:P1_S1]                           # set values of onshore  wind, season 1
    PR[1,2,:] = df_oper_params_R[:,:P1_S2]                           # set values of onshore  wind, season 2
    PR[1,3,:] = df_oper_params_R[:,:P1_S3]                           # set values of onshore  wind, season 3
    PR[1,4,:] = df_oper_params_R[:,:P1_S4]                           # set values of onshore  wind, season 4
    PR[2,1,:] = df_oper_params_R[:,:P2_S1]                           # set values of offshore wind, season 1
    PR[2,2,:] = df_oper_params_R[:,:P2_S2]                           # set values of offshore wind, season 2
    PR[2,3,:] = df_oper_params_R[:,:P2_S3]                           # set values of offshore wind, season 3
    PR[2,4,:] = df_oper_params_R[:,:P2_S4]                           # set values of offshore wind, season 4
    PR[3,1,:] = df_oper_params_R[:,:P3_S1]                           # set values of pv solar     , season 1
    PR[3,2,:] = df_oper_params_R[:,:P3_S2]                           # set values of pv solar     , season 2
    PR[3,3,:] = df_oper_params_R[:,:P3_S3]                           # set values of pv solar     , season 3
    PR[3,4,:] = df_oper_params_R[:,:P3_S4]                           # set values of pv solar     , season 4

    df_oper_params_D = CSV.read("$(pwd())/data/df_oper_params_D.csv" ;delim=",",use_mmap=false)
    PD = zeros(S[end],Hs[end])                                       # create empty array of energy demand (MW)
    PD[1,:] = df_oper_params_D[:,:P_S1]                              # set values of energy demand, season 1
    PD[2,:] = df_oper_params_D[:,:P_S2]                              # set values of energy demand, season 2
    PD[3,:] = df_oper_params_D[:,:P_S3]                              # set values of energy demand, season 3
    PD[4,:] = df_oper_params_D[:,:P_S4]                              # set values of energy demand, season 4

    df_oper_params_O = CSV.read("$(pwd())/data/df_oper_params_O.csv" ;delim=",",use_mmap=false)
    νD      = df_oper_params_O[1,:νD     ]                           # scaling demand factor
    c_shed  = df_oper_params_O[1,:c_shed ]                           # shedding demand cost (£/MW)
    α       = df_oper_params_O[1,:α      ]                           # seasonal weight 
    c_co2   = df_oper_params_O[1,:c_co2  ]                           # CO₂ cost (£/tCO₂)
    co2_lim = df_oper_params_O[1,:co2_lim]                           # CO₂ yearly budget (tCO₂)
    ρco2    = df_oper_params_O[1,:ρco2   ]                           # scaling factor of CO₂ budget

    ms  = ms_type(P,1:NI0,1:NI)
    mp  = mp_type(κ,x_hist,x_max,c_inv,c_fix,πi0,πi,itoi0)
    ps  = ps_type(P,G,B,R,S,Hs,G0,gn,b0,r0)
    pp  = pp_type(c_OMvarG,c_fuelG,em_co2G,ηG,rampG,ηB,PB,PR,PD,νD,c_shed,c_co2,co2_lim,ρco2,α)
    unc = u_type(NI,P[end]+size(H)[2],P[end],size(H)[2],size(C)[2],H,C)

    return ms,mp,ps,pp,unc
end