module WeakSeedGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.geometry.parameters"
const GRP_NAME = "weak_seed"

const PDATA = get_eb_parameters()

"""
    WeakSeed <: AbstractParameterGroup

Parameter group for weak seed geometry parameters.

# Fields
- `w_seed::`[`ParameterFloat`](@ref): $(PDATA.w_seed.description)
- `x_seed::`[`ParameterFloat`](@ref): $(PDATA.x_seed.description)
- `y_seed::`[`ParameterFloat`](@ref): $(PDATA.y_seed.description)
- `iuse_weak_seed::`[`ParameterInt`](@ref): $(PDATA.iuse_weak_seed.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of 
  parameter objects

# Nested Dot Access
- `w_seed = $(ROOT_NAME).$(GRP_NAME).w_seed.value`
- `x_seed = $(ROOT_NAME).$(GRP_NAME).x_seed.value`
- `y_seed = $(ROOT_NAME).$(GRP_NAME).y_seed.value`
- `iuse_weak_seed = $(ROOT_NAME).$(GRP_NAME).iuse_weak_seed.value`

# Constructor
    WeakSeed()

Create a new WeakSeed parameter group with default values.

# Returns
- `WeakSeed`: New WeakSeed parameter group with initialized default values

"""
mutable struct WeakSeed <: AbstractParameterGroup
    w_seed::ParameterFloat
    x_seed::ParameterFloat
    y_seed::ParameterFloat
    iuse_weak_seed::ParameterInt
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function WeakSeed()::WeakSeed
    pdata = get_eb_parameters()
    data = WeakSeed(
        pdata.w_seed,
        pdata.x_seed,
        pdata.y_seed,
        pdata.iuse_weak_seed,
        Union{ParameterFloat, ParameterInt}[]
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module 