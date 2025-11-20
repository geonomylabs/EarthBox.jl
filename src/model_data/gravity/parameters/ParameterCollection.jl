module ParameterCollection

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterCollection

const ROOT_NAME = "model.gravity.parameters"

const PDATA = get_eb_parameters()

"""
    Parameters <: AbstractParameterCollection

Parameter collection for gravity parameters.

# Fields
- `gravity_x::`[`ParameterFloat`](@ref): $(PDATA.gravity_x.description)
- `gravity_y::`[`ParameterFloat`](@ref): $(PDATA.gravity_y.description)
- `turn_off_gravity_y::`[`ParameterInt`](@ref): $(PDATA.turn_off_gravity_y.description)
- `nsteps_turn_off_gravity::`[`ParameterInt`](@ref): $(PDATA.nsteps_turn_off_gravity.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of parameter objects

# Nested Dot Access
- `gravity_x = $(ROOT_NAME).gravity_x.value`
- `gravity_y = $(ROOT_NAME).gravity_y.value`
- `turn_off_gravity_y = $(ROOT_NAME).turn_off_gravity_y.value`
- `nsteps_turn_off_gravity = $(ROOT_NAME).nsteps_turn_off_gravity.value`

# Constructor
    Parameters()

# Returns
- `Parameters`: New Parameters collection with initialized values
"""
mutable struct Parameters <: AbstractParameterCollection
    gravity_x::ParameterFloat
    gravity_y::ParameterFloat
    turn_off_gravity_y::ParameterInt
    nsteps_turn_off_gravity::ParameterInt
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function Parameters()::Parameters
    pdata = get_eb_parameters()
    data = Parameters(
        pdata.gravity_x,
        pdata.gravity_y,
        pdata.turn_off_gravity_y,
        pdata.nsteps_turn_off_gravity,
        Union{ParameterFloat, ParameterInt}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
