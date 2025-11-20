module CrustalHoleGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.geometry.parameters"
const GRP_NAME = "crustal_hole"

const PDATA = get_eb_parameters()

"""
    CrustalHole <: AbstractParameterGroup

Parameter group for crustal hole geometry parameters.

# Fields
- `xhole_start::`[`ParameterFloat`](@ref): $(PDATA.xhole_start.description)
- `xhole_middle::`[`ParameterFloat`](@ref): $(PDATA.xhole_middle.description)
- `xhole_end::`[`ParameterFloat`](@ref): $(PDATA.xhole_end.description)
- `xhole_depth::`[`ParameterFloat`](@ref): $(PDATA.xhole_depth.description)
- `obj_list::Vector{Union{[`ParameterFloat`](@ref)}}`: List of parameter objects

# Nested Dot Access
- `xhole_start = $(ROOT_NAME).$(GRP_NAME).xhole_start.value`
- `xhole_middle = $(ROOT_NAME).$(GRP_NAME).xhole_middle.value`
- `xhole_end = $(ROOT_NAME).$(GRP_NAME).xhole_end.value`
- `xhole_depth = $(ROOT_NAME).$(GRP_NAME).xhole_depth.value`

# Constructor
    CrustalHole()

Create a new CrustalHole parameter group with default values.

# Returns
- `CrustalHole`: New CrustalHole parameter group with initialized default values

"""
mutable struct CrustalHole <: AbstractParameterGroup
    xhole_start::ParameterFloat
    xhole_middle::ParameterFloat
    xhole_end::ParameterFloat
    xhole_depth::ParameterFloat
    obj_list::Vector{Union{ParameterFloat}}
end

function CrustalHole()::CrustalHole
    pdata = get_eb_parameters()
    data = CrustalHole(
        pdata.xhole_start,
        pdata.xhole_middle,
        pdata.xhole_end,
        pdata.xhole_depth,
        Union{ParameterFloat}[]
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
