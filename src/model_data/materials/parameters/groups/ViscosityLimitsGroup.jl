module ViscosityLimitsGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.materials.parameters"
const GRP_NAME = "viscosity_limits"

const PDATA = get_eb_parameters()

"""
    ViscosityLimits <: AbstractParameterGroup

Parameter group for viscosity limits.

# Fields
- `viscosity_min::`[`ParameterFloat`](@ref): $(PDATA.viscosity_min.description)
- `viscosity_max::`[`ParameterFloat`](@ref): $(PDATA.viscosity_max.description)
- `obj_list::Vector{ParameterFloat}`: List of parameter objects

# Nested Dot Access
- `viscosity_min = $(ROOT_NAME).$(GRP_NAME).viscosity_min.value`
- `viscosity_max = $(ROOT_NAME).$(GRP_NAME).viscosity_max.value`

# Constructor
    ViscosityLimits()

# Returns
- `ViscosityLimits`: New ViscosityLimits parameter group with initialized values

"""
mutable struct ViscosityLimits <: AbstractParameterGroup
    viscosity_min::ParameterFloat
    viscosity_max::ParameterFloat
    obj_list::Vector{ParameterFloat}
end

function ViscosityLimits()::ViscosityLimits
    pdata = get_eb_parameters()
    data = ViscosityLimits(
        pdata.viscosity_min,
        pdata.viscosity_max,
        ParameterFloat[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
