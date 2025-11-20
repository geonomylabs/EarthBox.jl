module RandomFrictionGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.materials.parameters"
const GRP_NAME = "random_friction"

const PDATA = get_eb_parameters()

"""
    RandomFriction <: AbstractParameterGroup

Parameter group for random friction properties.

# Fields
- `iuse_random_fric::`[`ParameterInt`](@ref): $(PDATA.iuse_random_fric.description)
- `delta_fric_coef::`[`ParameterFloat`](@ref): $(PDATA.delta_fric_coef.description)
- `iuse_random_fric_time::`[`ParameterInt`](@ref): $(PDATA.iuse_random_fric_time.description)
- `randomization_factor::`[`ParameterFloat`](@ref): $(PDATA.randomization_factor.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of parameter objects

# Nested Dot Access
- `iuse_random_fric = $(ROOT_NAME).$(GRP_NAME).iuse_random_fric.value`
- `delta_fric_coef = $(ROOT_NAME).$(GRP_NAME).delta_fric_coef.value`
- `iuse_random_fric_time = $(ROOT_NAME).$(GRP_NAME).iuse_random_fric_time.value`
- `randomization_factor = $(ROOT_NAME).$(GRP_NAME).randomization_factor.value`

# Constructor
    RandomFriction()

Create a new RandomFriction parameter group with default values.

# Returns
- `RandomFriction`: New RandomFriction parameter group with initialized values

"""
mutable struct RandomFriction <: AbstractParameterGroup
    iuse_random_fric::ParameterInt
    delta_fric_coef::ParameterFloat
    iuse_random_fric_time::ParameterInt
    randomization_factor::ParameterFloat
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function RandomFriction()::RandomFriction
    pdata = get_eb_parameters()
    data = RandomFriction(
        pdata.iuse_random_fric,
        pdata.delta_fric_coef,
        pdata.iuse_random_fric_time,
        pdata.randomization_factor,
        Union{ParameterFloat, ParameterInt}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
