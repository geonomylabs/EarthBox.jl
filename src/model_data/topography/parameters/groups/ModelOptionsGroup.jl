module ModelOptionsGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt, ParameterStr
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.topography.parameters"
const GRP_NAME = "model_options"

const PDATA = get_eb_parameters()

"""
    ModelOptions <: AbstractParameterGroup

Parameter group for topography model options.

# Fields
- `iuse_topo::`[`ParameterInt`](@ref): $(PDATA.iuse_topo.description)
- `itype_topo::`[`ParameterInt`](@ref): $(PDATA.itype_topo.description)
- `stype_topo::`[`ParameterStr`](@ref): $(PDATA.stype_topo.description)
- `iuse_downhill_diffusion::`[`ParameterInt`](@ref): $(PDATA.iuse_downhill_diffusion.description)
- `iuse_salt_deposition::`[`ParameterInt`](@ref): $(PDATA.iuse_salt_deposition.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt, ParameterStr}}`: List of parameter objects

# Nested Dot Access
- `iuse_topo = $(ROOT_NAME).$(GRP_NAME).iuse_topo.value`
- `itype_topo = $(ROOT_NAME).$(GRP_NAME).itype_topo.value`
- `stype_topo = $(ROOT_NAME).$(GRP_NAME).stype_topo.value`
- `iuse_downhill_diffusion = $(ROOT_NAME).$(GRP_NAME).iuse_downhill_diffusion.value`
- `iuse_salt_deposition = $(ROOT_NAME).$(GRP_NAME).iuse_salt_deposition.value`

# Constructor
    ModelOptions()

# Returns
- `ModelOptions`: New ModelOptions parameter group with initialized values

"""
mutable struct ModelOptions <: AbstractParameterGroup
    iuse_topo::ParameterInt
    itype_topo::ParameterInt
    stype_topo::ParameterStr
    iuse_downhill_diffusion::ParameterInt
    iuse_salt_deposition::ParameterInt
    obj_list::Vector{Union{ParameterFloat, ParameterInt, ParameterStr}}
end

function ModelOptions()::ModelOptions
    pdata = get_eb_parameters()
    data = ModelOptions(
        pdata.iuse_topo,
        pdata.itype_topo,
        pdata.stype_topo,
        pdata.iuse_downhill_diffusion,
        pdata.iuse_salt_deposition,
        Union{ParameterFloat, ParameterInt, ParameterStr}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
