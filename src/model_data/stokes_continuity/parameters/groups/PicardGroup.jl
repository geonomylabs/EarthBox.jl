module PicardGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt, ParameterStr
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.stokes_continuity.parameters"
const GRP_NAME = "picard"

const PDATA = get_eb_parameters()

"""
    Picard <: AbstractParameterGroup

Parameter group for Picard iteration parameters.

# Fields
- `nglobal::`[`ParameterInt`](@ref): $(PDATA.nglobal.description)
- `tolerance_picard::`[`ParameterFloat`](@ref): $(PDATA.tolerance_picard.description)
- `iconverge::`[`ParameterInt`](@ref): $(PDATA.iconverge.description)
- `iglobal::`[`ParameterInt`](@ref): $(PDATA.iglobal.description)
- `itype_global::`[`ParameterInt`](@ref): $(PDATA.itype_global.description)
- `stype_global::`[`ParameterStr`](@ref): $(PDATA.stype_global.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of parameter objects

# Nested Dot Access
- `nglobal = $(ROOT_NAME).$(GRP_NAME).nglobal.value`
- `tolerance_picard = $(ROOT_NAME).$(GRP_NAME).tolerance_picard.value`
- `iconverge = $(ROOT_NAME).$(GRP_NAME).iconverge.value`
- `iglobal = $(ROOT_NAME).$(GRP_NAME).iglobal.value`
- `itype_global = $(ROOT_NAME).$(GRP_NAME).itype_global.value`
- `stype_global = $(ROOT_NAME).$(GRP_NAME).stype_global.value`

# Constructor
    Picard()

Create a new Picard parameter group with default values.

# Returns
- `Picard`: New Picard parameter group with initialized values

"""
mutable struct Picard <: AbstractParameterGroup
    nglobal::ParameterInt
    tolerance_picard::ParameterFloat
    iconverge::ParameterInt
    iglobal::ParameterInt
    itype_global::ParameterInt
    stype_global::ParameterStr
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function Picard()::Picard
    pdata = get_eb_parameters()
    data = Picard(
        pdata.nglobal,
        pdata.tolerance_picard,
        pdata.iconverge,
        pdata.iglobal,
        pdata.itype_global,
        pdata.stype_global,
        Union{ParameterFloat, ParameterInt}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
