module MeltingMaterialIDsGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterInt
import EarthBox.Parameters: print_parameter
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.melting.parameters"
const GRP_NAME = "melting_material_ids"

const PDATA = get_eb_parameters()

"""
    MeltingMaterialIDs <: AbstractParameterGroup

Parameter group for material IDs used in melting processes.

# Fields
- `ipmf1::`[`ParameterInt`](@ref): $(PDATA.ipmf1.description)
- `ipmf2::`[`ParameterInt`](@ref): $(PDATA.ipmf2.description)
- `ipmf3::`[`ParameterInt`](@ref): $(PDATA.ipmf3.description)
- `ipmf4::`[`ParameterInt`](@ref): $(PDATA.ipmf4.description)
- `ipmf_molten::`[`ParameterInt`](@ref): $(PDATA.ipmf_molten.description)
- `ipmf_solid::`[`ParameterInt`](@ref): $(PDATA.ipmf_solid.description)
- `ipmf_molten_vol::`[`ParameterInt`](@ref): $(PDATA.ipmf_molten_vol.description)
- `ipmf_solid_vol::`[`ParameterInt`](@ref): $(PDATA.ipmf_solid_vol.description)
- `obj_list::Vector{ParameterInt}`: List of parameter objects

# Nested Dot Access
- `ipmf1 = $(ROOT_NAME).$(GRP_NAME).ipmf1.value`
- `ipmf2 = $(ROOT_NAME).$(GRP_NAME).ipmf2.value`
- `ipmf3 = $(ROOT_NAME).$(GRP_NAME).ipmf3.value`
- `ipmf4 = $(ROOT_NAME).$(GRP_NAME).ipmf4.value`
- `ipmf_molten = $(ROOT_NAME).$(GRP_NAME).ipmf_molten.value`
- `ipmf_solid = $(ROOT_NAME).$(GRP_NAME).ipmf_solid.value`
- `ipmf_molten_vol = $(ROOT_NAME).$(GRP_NAME).ipmf_molten_vol.value`
- `ipmf_solid_vol = $(ROOT_NAME).$(GRP_NAME).ipmf_solid_vol.value`

# Constructor
    MeltingMaterialIDs()

# Returns
- `MeltingMaterialIDs`: New MeltingMaterialIDs parameter group with initialized values

"""
mutable struct MeltingMaterialIDs <: AbstractParameterGroup
    ipmf1::ParameterInt
    ipmf2::ParameterInt
    ipmf3::ParameterInt
    ipmf4::ParameterInt
    ipmf_molten::ParameterInt
    ipmf_solid::ParameterInt
    ipmf_molten_vol::ParameterInt
    ipmf_solid_vol::ParameterInt
    obj_list::Vector{ParameterInt}
end

function MeltingMaterialIDs()::MeltingMaterialIDs
    pdata = get_eb_parameters()
    data = MeltingMaterialIDs(
        pdata.ipmf1,
        pdata.ipmf2,
        pdata.ipmf3,
        pdata.ipmf4,
        pdata.ipmf_molten,
        pdata.ipmf_solid,
        pdata.ipmf_molten_vol,
        pdata.ipmf_solid_vol,
        ParameterInt[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
