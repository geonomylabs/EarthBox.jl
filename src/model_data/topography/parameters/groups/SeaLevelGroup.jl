module SeaLevelGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt, ParameterStr
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.topography.parameters"
const GRP_NAME = "sealevel"

const PDATA = get_eb_parameters()

"""
    SeaLevel <: AbstractParameterGroup

Parameter group for sea level parameters.

# Fields
- `itype_sealevel::`[`ParameterInt`](@ref): $(PDATA.itype_sealevel.description)
- `stype_sealevel::`[`ParameterStr`](@ref): $(PDATA.stype_sealevel.description)
- `y_water_ini::`[`ParameterFloat`](@ref): $(PDATA.y_water_ini.description)
- `base_level_shift::`[`ParameterFloat`](@ref): $(PDATA.base_level_shift.description)
- `base_level_shift_end_time::`[`ParameterFloat`](@ref): $(PDATA.base_level_shift_end_time.description)
- `y_sealevel::`[`ParameterFloat`](@ref): $(PDATA.y_sealevel.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt, ParameterStr}}`: List of parameter objects

# Nested Dot Access
- `itype_sealevel = $(ROOT_NAME).$(GRP_NAME).itype_sealevel.value`
- `stype_sealevel = $(ROOT_NAME).$(GRP_NAME).stype_sealevel.value`
- `y_water_ini = $(ROOT_NAME).$(GRP_NAME).y_water_ini.value`
- `base_level_shift = $(ROOT_NAME).$(GRP_NAME).base_level_shift.value`
- `base_level_shift_end_time = $(ROOT_NAME).$(GRP_NAME).base_level_shift_end_time.value`
- `y_sealevel = $(ROOT_NAME).$(GRP_NAME).y_sealevel.value`

# Constructor
    SeaLevel()

# Returns
- `SeaLevel`: New SeaLevel parameter group with initialized values

"""
mutable struct SeaLevel <: AbstractParameterGroup
    itype_sealevel::ParameterInt
    stype_sealevel::ParameterStr
    y_water_ini::ParameterFloat
    base_level_shift::ParameterFloat
    base_level_shift_end_time::ParameterFloat
    y_sealevel::ParameterFloat
    obj_list::Vector{Union{ParameterFloat, ParameterInt, ParameterStr}}
end

function SeaLevel()::SeaLevel
    pdata = get_eb_parameters()
    data = SeaLevel(
        pdata.itype_sealevel,
        pdata.stype_sealevel,
        pdata.y_water_ini,
        pdata.base_level_shift,
        pdata.base_level_shift_end_time,
        pdata.y_sealevel,
        Union{ParameterFloat, ParameterInt, ParameterStr}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
