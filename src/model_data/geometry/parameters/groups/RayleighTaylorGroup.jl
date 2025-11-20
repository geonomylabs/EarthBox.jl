module RayleighTaylorGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.geometry.parameters"
const GRP_NAME = "rayleigh_taylor"

const PDATA = get_eb_parameters()

"""
    RayleighTaylor <: AbstractParameterGroup

Parameter group for Rayleigh-Taylor instability geometry parameters.

# Fields
- `depth_interface_h1::`[`ParameterFloat`](@ref): $(PDATA.depth_interface_h1.description)
- `wave_length_lambda::`[`ParameterFloat`](@ref): $(PDATA.wave_length_lambda.description)
- `amplitude_initial::`[`ParameterFloat`](@ref): $(PDATA.amplitude_initial.description)
- `obj_list::Vector{ParameterFloat}`: List of parameter objects

# Nested Dot Access
- `depth_interface_h1 = $(ROOT_NAME).$(GRP_NAME).depth_interface_h1.value`
- `wave_length_lambda = $(ROOT_NAME).$(GRP_NAME).wave_length_lambda.value`
- `amplitude_initial = $(ROOT_NAME).$(GRP_NAME).amplitude_initial.value`

# Constructor
    RayleighTaylor()

Create a new RayleighTaylor parameter group with default values.

# Returns
- `RayleighTaylor`: New RayleighTaylor parameter group with initialized default values

"""
mutable struct RayleighTaylor <: AbstractParameterGroup
    depth_interface_h1::ParameterFloat
    wave_length_lambda::ParameterFloat
    amplitude_initial::ParameterFloat
    obj_list::Vector{ParameterFloat}
end

function RayleighTaylor()::RayleighTaylor
    pdata = get_eb_parameters()
    data = RayleighTaylor(
        pdata.depth_interface_h1,
        pdata.wave_length_lambda,
        pdata.amplitude_initial,
        ParameterFloat[]
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
