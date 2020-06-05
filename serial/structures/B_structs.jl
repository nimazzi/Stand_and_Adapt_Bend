module B_structs

using JuMP
using O_structs

# */ --------------------------------------------------------------------------------------------- /* #
# */ --- structures of data : on main core ------------------------------------------------------- /* #
# */ --------------------------------------------------------------------------------------------- /* #

export D__type
export R__type
export T1_type
export T2_type
export T3_type
export H__type
export M__type
export B1_type
export B2_type
export B3_type

mutable struct D__type
    ms::ms_type
    mp::mp_type
    ps::ps_type
    pp::pp_type
   unc::u_type
  case::Int64
  algm::Int64
     q::Float64
     w::Int64
  kidx::Array{BitArray{1},1}
  iidx::Array{Array{Int64,1},1}
     J::Int64
     δ::Float64
end

mutable struct R__type
     m::JuMP.Model
  vars::Array{JuMP.VariableRef,2}
  cons::Array{JuMP.ConstraintRef,2}
end

mutable struct T1_type
     x::Array{Float64,2}
     c::Array{Float64,2}
     β::Array{Float64,1}
     θ::Array{Float64,1}
     λ::Array{Float64,2}
     Δ::Float64
end

mutable struct T2_type
     x::Array{Float64,2}
     c::Array{Float64,2}
     β::Array{Float64,1}
    θl::Array{Float64,1}
    λl::Array{Float64,2}
    θu::Array{Float64,1}
    ϕu::Array{Float64,2}
    Ie::Array{Int64,1}
     Δ::Float64
end

mutable struct T3_type
     x::Array{Float64,2}
     c::Array{Float64,2}
     β::Array{Float64,1}
    θp::Array{Float64,1}
    θd::Array{Float64,1}
     λ::Array{Float64,2}
     Δ::Float64
end

mutable struct H__type
     L::Array{Float64,1}
     U::Array{Float64,1}
     t::Array{Float64,1}
     T::Array{Float64,2}
     k::Int64
end

mutable struct M__type
      n::Int64
      θ::Array{Float64,1}
      λ::Array{Float64,2}
      ϕ::Array{Float64,2}
      x::Array{Float64,2}
      c::Array{Float64,2}
end

mutable struct B1_type
   rmp::R__type
  temp::T1_type
  hist::H__type
  data::D__type
end

mutable struct B2_type
     rmp::R__type
    temp::T2_type
    hist::H__type
      m0::M__type
    data::D__type
 end

mutable struct B3_type
     rmp::R__type
    temp::T3_type
    hist::H__type
    data::D__type
end

# */ --------------------------------------------------------------------------------------------- /* #
# */ --- structures of data : everywhere --------------------------------------------------------- /* #
# */ --------------------------------------------------------------------------------------------- /* #

export E__type
export O__type
export d__type
export t1_type
export t2_type
export t3_type
export S1_type
export S2_type
export S3_type

mutable struct E__type
      m::JuMP.Model
   vars::Array{JuMP.VariableRef,1}
   objf::GenericAffExpr{Float64,VariableRef}
end

mutable struct O__type
      m::JuMP.Model
   cons::Array{JuMP.ConstraintRef,1}
   vars::Array{JuMP.VariableRef,1}
   objf::GenericAffExpr{Float64,VariableRef}
   cexp::GenericAffExpr{Float64,VariableRef}
   help::Array{Float64,1}
      n::Int64
end

mutable struct d__type
     ps::ps_type
     pp::pp_type
    unc::u_type
   algm::Int64
      J::Int64
end

mutable struct t1_type
      x::Array{Float64,2}
      c::Array{Float64,2}
      θ::Array{Float64,1}
      λ::Array{Float64,2}
      i::Int64
end

mutable struct t2_type
     x::Array{Float64,2}
     c::Array{Float64,2}
    θl::Array{Float64,1}
    θu::Array{Float64,1}
    λl::Array{Float64,2}
    ϕu::Array{Float64,2}
     i::Int64
end

mutable struct t3_type
     x::Array{Float64,2}
     c::Array{Float64,2}
    θp::Array{Float64,1}
    θd::Array{Float64,1}
     λ::Array{Float64,2}
     i::Int64
end

mutable struct S1_type
     ex::E__type
   data::d__type
   temp::t1_type
end

mutable struct S2_type
     ex::E__type
    olb::O__type
    oub::O__type
      m::M__type
   data::d__type
   temp::t2_type
end

mutable struct S3_type
     ex::E__type
   data::d__type
   temp::t3_type
end

# */ --------------------------------------------------------------------------------------------- /* #
# */ --------------------------------------------------------------------------------------------- /* #

end
