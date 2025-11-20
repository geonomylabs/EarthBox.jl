module MaterialDescriptionGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterInt, ParameterStr
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.materials.parameters"
const GRP_NAME = "description"

const PDATA = get_eb_parameters()

"""
    MaterialDescription <: AbstractParameterGroup

Parameter group for material and plasticity model descriptions.

# Fields
- `itype_mat::`[`ParameterInt`](@ref): $(PDATA.itype_mat.description)
- `stype_mat::`[`ParameterStr`](@ref): $(PDATA.stype_mat.description)
- `nmats::`[`ParameterInt`](@ref): $(PDATA.nmats.description)
- `itype_plasticity::`[`ParameterInt`](@ref): $(PDATA.itype_plasticity.description)
- `stype_plasticity::`[`ParameterStr`](@ref): $(PDATA.stype_plasticity.description)
- `obj_list::Vector{ParameterInt}`: List of parameter objects

# Nested Dot Access
- `itype_mat = $(ROOT_NAME).$(GRP_NAME).itype_mat.value`
- `stype_mat = $(ROOT_NAME).$(GRP_NAME).stype_mat.value`
- `nmats = $(ROOT_NAME).$(GRP_NAME).nmats.value`
- `itype_plasticity = $(ROOT_NAME).$(GRP_NAME).itype_plasticity.value`
- `stype_plasticity = $(ROOT_NAME).$(GRP_NAME).stype_plasticity.value`

# Constructor
    MaterialDescription()

Create a new MaterialDescription parameter group with default values.

# Returns
- `MaterialDescription`: New MaterialDescription parameter group with initialized values

"""
mutable struct MaterialDescription <: AbstractParameterGroup
    itype_mat::ParameterInt
    stype_mat::ParameterStr
    nmats::ParameterInt
    itype_plasticity::ParameterInt
    stype_plasticity::ParameterStr
    obj_list::Vector{ParameterInt}
end

function MaterialDescription()::MaterialDescription
    pdata = get_eb_parameters()
    data = MaterialDescription(
        pdata.itype_mat,
        pdata.stype_mat,
        pdata.nmats,
        pdata.itype_plasticity,
        pdata.stype_plasticity,
        ParameterInt[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
