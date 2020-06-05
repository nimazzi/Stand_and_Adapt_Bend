module O_structs

export u_type,ms_type,mp_type,ps_type,pp_type

mutable struct u_type
      ni::Int64                  # number of subproblems
      nx::Int64                  # size of vector x
     nx0::Int64                  # size of master vector decisions x0
      nh::Int64                  # size of vector h
      nc::Int64                  # size of vector c
       h::Array{Float64,2}       # array of vector h
       c::Array{Float64,2}       # array of vector c
end

mutable struct ms_type 
       P::UnitRange{Int64}       # set of technologies
      I0::UnitRange{Int64}       # set of "investment" nodes
       I::UnitRange{Int64}       # set of "operational" nodes
end

mutable struct mp_type 
        κ::Float64                # years of operational problem
       xh::Array{Float64,2}       # historical installed capacity (GW)
       xm::Array{Float64,1}       # maximum accumulated capacity (GW)
       ci::Array{Float64,2}       # investment cost (£/GW)
       cf::Array{Float64,1}       # fix OM cost (£/GWyr)
       π0::Array{Float64,1}       # probability associated to "investment" node i0 (-)
        π::Array{Float64,1}       # probability associated to "operational" node i (-)
      map::Vector{Array{Int64,1}} # map "operational" node i → "investment" node i0
end

mutable struct ps_type
        P::UnitRange{Int64}       # set of technologies
        G::UnitRange{Int64}       # set of conventional technologies
        B::UnitRange{Int64}       # set of storage technologies
        R::UnitRange{Int64}       # set of renewable technologies
        S::UnitRange{Int64}       # set of seasons
        H::UnitRange{Int64}       # hours within a season
       G0::UnitRange{Int64}       # set of conventional technologies exept nuclear
       gn::Int64                  # nuclear index
       b0::Int64                  # index 0 in x of storage technologies
       r0::Int64                  # index 0 in x of renewables technologies
end

mutable struct pp_type
      cvg::Array{Float64,1}       # var OM cost (conv. tech.) (£/MWh)
      cfg::Array{Float64,1}       # fuel cost (conv. tech.) (£/MWh)
       eg::Array{Float64,1}       # CO2 emissions (conv. tech.) (tCO₂/MWh)
       ηg::Array{Float64,1}       # efficiency (conv. tech.) (-)
       rg::Array{Float64,1}       # ramping limitations (conv. tech.) (MW)
       ηb::Array{Float64,1}       # efficiency (storage tech.) (-)
       pb::Array{Float64,1}       # ch/disch power (storage tech.) (MW/MW)
       pr::Array{Float64,3}       # renewable energy prod (re. tech.) (MW/MW)
       pd::Array{Float64,2}       # energy demand (MW)
       νd::Float64                # scaling demand (-)
       cs::Float64                # shedding demand cost (£/MW)
       co::Float64                # CO2 cost (conv. tech.) (£/tCO₂)
      lco::Float64                # co2 yearly limit (tCO₂)
      ρco::Float64                # scaling co2 limit (-)
        α::Float64                # seasonal weight (-)
end

end