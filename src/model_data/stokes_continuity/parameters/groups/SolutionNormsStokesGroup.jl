module SolutionNormsStokesGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.stokes_continuity.parameters"
const GRP_NAME = "solution_norms"

const PDATA = get_eb_parameters()

"""
    SolutionNormsStokes <: AbstractParameterGroup

Parameter group for Stokes solution norm tracking.

# Fields
- `dsoluv1_abs_inf::`[`ParameterFloat`](@ref): $(PDATA.dsoluv1_abs_inf.description)
- `dsoluv1_rel_inf::`[`ParameterFloat`](@ref): $(PDATA.dsoluv1_rel_inf.description)
- `dsoluv1_abs_L2::`[`ParameterFloat`](@ref): $(PDATA.dsoluv1_abs_L2.description)
- `dsoluv1_rel_L2::`[`ParameterFloat`](@ref): $(PDATA.dsoluv1_rel_L2.description)
- `dvx1_abs_L2::`[`ParameterFloat`](@ref): $(PDATA.dvx1_abs_L2.description)
- `dvx1_rel_L2::`[`ParameterFloat`](@ref): $(PDATA.dvx1_rel_L2.description)
- `dvy1_abs_L2::`[`ParameterFloat`](@ref): $(PDATA.dvy1_abs_L2.description)
- `dvy1_rel_L2::`[`ParameterFloat`](@ref): $(PDATA.dvy1_rel_L2.description)
- `dpr1_abs_L2::`[`ParameterFloat`](@ref): $(PDATA.dpr1_abs_L2.description)
- `dpr1_rel_L2::`[`ParameterFloat`](@ref): $(PDATA.dpr1_rel_L2.description)
- `dvxy_abs_inf::`[`ParameterFloat`](@ref): $(PDATA.dvxy_abs_inf.description)
- `dvxy_rel_inf::`[`ParameterFloat`](@ref): $(PDATA.dvxy_rel_inf.description)
- `dvxy_abs_L2::`[`ParameterFloat`](@ref): $(PDATA.dvxy_abs_L2.description)
- `dvxy_rel_L2::`[`ParameterFloat`](@ref): $(PDATA.dvxy_rel_L2.description)
- `global_yield_error::`[`ParameterFloat`](@ref): $(PDATA.global_yield_error.description)
- `obj_list::Vector{ParameterFloat}`: List of parameter objects

# Nested Dot Access
- `dsoluv1_abs_inf = $(ROOT_NAME).$(GRP_NAME).dsoluv1_abs_inf.value`
- `dsoluv1_abs_L2 = $(ROOT_NAME).$(GRP_NAME).dsoluv1_abs_L2.value`
- `dvx1_abs_L2 = $(ROOT_NAME).$(GRP_NAME).dvx1_abs_L2.value`

# Constructor
    SolutionNormsStokes()

# Returns
- `SolutionNormsStokes`: New SolutionNormsStokes parameter group with initialized values

"""
mutable struct SolutionNormsStokes <: AbstractParameterGroup
    dsoluv1_abs_inf::ParameterFloat
    dsoluv1_rel_inf::ParameterFloat
    dsoluv1_abs_L2::ParameterFloat
    dsoluv1_rel_L2::ParameterFloat
    dvx1_abs_L2::ParameterFloat
    dvx1_rel_L2::ParameterFloat
    dvy1_abs_L2::ParameterFloat
    dvy1_rel_L2::ParameterFloat
    dpr1_abs_L2::ParameterFloat
    dpr1_rel_L2::ParameterFloat
    dvxy_abs_inf::ParameterFloat
    dvxy_rel_inf::ParameterFloat
    dvxy_abs_L2::ParameterFloat
    dvxy_rel_L2::ParameterFloat
    global_yield_error::ParameterFloat
    obj_list::Vector{ParameterFloat}
end

function SolutionNormsStokes()::SolutionNormsStokes
    pdata = get_eb_parameters()
    data = SolutionNormsStokes(
        pdata.dsoluv1_abs_inf,
        pdata.dsoluv1_rel_inf,
        pdata.dsoluv1_abs_L2,
        pdata.dsoluv1_rel_L2,
        pdata.dvx1_abs_L2,
        pdata.dvx1_rel_L2,
        pdata.dvy1_abs_L2,
        pdata.dvy1_rel_L2,
        pdata.dpr1_abs_L2,
        pdata.dpr1_rel_L2,
        pdata.dvxy_abs_inf,
        pdata.dvxy_rel_inf,
        pdata.dvxy_abs_L2,
        pdata.dvxy_rel_L2,
        pdata.global_yield_error,
        ParameterFloat[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
