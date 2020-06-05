function load_data(case::Int64)::Tuple{ms_type,mp_type,ps_type,pp_type,u_type}

    df_sets = CSV.read("$(pwd())/data/df_sets.csv" ;delim=",",use_mmap=false)
    P  = 1:df_sets[1,:P]
    G  = 1:df_sets[1,:G]
    B  = 1:df_sets[1,:B]
    R  = 1:df_sets[1,:R]
    S  = 1:df_sets[1,:S]
    Hs = 1:df_sets[1,:H]
    G0 = 1:df_sets[1,:G0]
    gn = df_sets[1,:gn]
    b0 = df_sets[1,:b0]
    r0 = df_sets[1,:r0]

    df_unc_sets = CSV.read("$(pwd())/data/df_unc_sets_case$(case).csv" ;delim=",",use_mmap=false)
    NI1 = df_unc_sets[1,:NΩ1]
    NI2 = df_unc_sets[1,:NΩ2]
    NI3 = df_unc_sets[1,:NΩ3]
    NI0 = NI1 + NI2
    NI  = NI2 + NI3

    df_unc_params = CSV.read("$(pwd())/data/df_unc_params_case$(case).csv" ;delim=",",use_mmap=false)
    H      = zeros(NI,2)
    C      = zeros(NI,2)
    H[:,1] = df_unc_params[:,:ρco2 ]
    H[:,2] = -df_unc_params[:,:νD   ]
    C[:,1] = df_unc_params[:,:c_co2]
    C[:,2] = df_unc_params[:,:c_ur ]

    df_inv_params_P = CSV.read("$(pwd())/data/df_inv_params_P.csv" ;delim=",",use_mmap=false)
    c_inv_0   = df_inv_params_P[:,:c_inv_0]
    c_inv_5   = df_inv_params_P[:,:c_inv_5  ]
    c_fix     = df_inv_params_P[:,:c_fix    ].*10^3
    x_hist_5  = df_inv_params_P[:,:x_hist_5 ]
    x_hist_10 = df_inv_params_P[:,:x_hist_10]
    x_max     = df_inv_params_P[:,:x_max    ].*10^-3
    c_inv     = hcat(repeat(c_inv_0 ,1,NI1),repeat(c_inv_5  ,1,NI2)).*10^3
    x_hist    = hcat(repeat(x_hist_5,1,NI2),repeat(x_hist_10,1,NI3)).*10^-3
    πi0       = vcat(1/NI1*ones(NI1),1/NI2*ones(NI2))
    πi        = vcat(1/NI2*ones(NI2),1/NI3*ones(NI3))
    itoi0     = Vector{Array{Int64,1}}()  # [ωM,..] = ΩStoM[ωS]
    for i in 1:NI2 push!(itoi0,[1]) end
    for i in (NI2+1):(NI2+NI3) push!(itoi0,[1,Int64(((i-NI2-1)-(i-NI2-1)%NI2)/NI2+1)+1]) end

    df_inv_params_O = CSV.read("$(pwd())/data/df_inv_params_O.csv" ;delim=",",use_mmap=false)
    κ = df_inv_params_O[1,:κ]

    df_oper_params_G = CSV.read("$(pwd())/data/df_oper_params_G.csv" ;delim=",",use_mmap=false)
    c_OMvarG = df_oper_params_G[:,:c_varOM]
    c_fuelG  = df_oper_params_G[:,:c_fuel ]
    em_co2G  = df_oper_params_G[:,:em_co2 ]
    ηG       = df_oper_params_G[:,:η      ]
    rampG    = df_oper_params_G[:,:ramp   ]

    df_oper_params_B = CSV.read("$(pwd())/data/df_oper_params_B.csv" ;delim=",",use_mmap=false)
    ηB = df_oper_params_B[:,:η]
    PB = df_oper_params_B[:,:P]

    df_oper_params_R = CSV.read("$(pwd())/data/df_oper_params_R.csv" ;delim=",",use_mmap=false)
    PR = zeros(R[end],S[end],Hs[end])
    PR[1,1,:] = df_oper_params_R[:,:P1_S1]
    PR[1,2,:] = df_oper_params_R[:,:P1_S2]
    PR[1,3,:] = df_oper_params_R[:,:P1_S3]
    PR[1,4,:] = df_oper_params_R[:,:P1_S4]
    PR[2,1,:] = df_oper_params_R[:,:P2_S1]
    PR[2,2,:] = df_oper_params_R[:,:P2_S2]
    PR[2,3,:] = df_oper_params_R[:,:P2_S3]
    PR[2,4,:] = df_oper_params_R[:,:P2_S4]
    PR[3,1,:] = df_oper_params_R[:,:P3_S1]
    PR[3,2,:] = df_oper_params_R[:,:P3_S2]
    PR[3,3,:] = df_oper_params_R[:,:P3_S3]
    PR[3,4,:] = df_oper_params_R[:,:P3_S4]

    df_oper_params_D = CSV.read("$(pwd())/data/df_oper_params_D.csv" ;delim=",",use_mmap=false)
    PD = zeros(S[end],Hs[end])
    PD[1,:] = df_oper_params_D[:,:P_S1]
    PD[2,:] = df_oper_params_D[:,:P_S2]
    PD[3,:] = df_oper_params_D[:,:P_S3]
    PD[4,:] = df_oper_params_D[:,:P_S4]

    df_oper_params_O = CSV.read("$(pwd())/data/df_oper_params_O.csv" ;delim=",",use_mmap=false)
    νD      = df_oper_params_O[1,:νD     ]
    c_shed  = df_oper_params_O[1,:c_shed ]
    α       = df_oper_params_O[1,:α      ]
    c_co2   = df_oper_params_O[1,:c_co2  ]
    co2_lim = df_oper_params_O[1,:co2_lim]
    ρco2    = df_oper_params_O[1,:ρco2   ]

    ms  = ms_type(P,1:NI0,1:NI)
    mp  = mp_type(κ,x_hist,x_max,c_inv,c_fix,πi0,πi,itoi0)
    ps  = ps_type(P,G,B,R,S,Hs,G0,gn,b0,r0)
    pp  = pp_type(c_OMvarG,c_fuelG,em_co2G,ηG,rampG,ηB,PB,PR,PD,νD,c_shed,c_co2,co2_lim,ρco2,α)
    unc = u_type(NI,P[end]+size(H)[2],P[end],size(H)[2],size(C)[2],H,C)

    return ms,mp,ps,pp,unc
end