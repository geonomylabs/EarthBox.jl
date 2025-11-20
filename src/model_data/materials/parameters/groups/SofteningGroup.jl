module SofteningGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.materials.parameters"
const GRP_NAME = "softening"

const PDATA = get_eb_parameters()

"""
    Softening <: AbstractParameterGroup

Parameter group for viscous strain softening properties.

# Fields
- `iuse_viscous_strain_soft::`[`ParameterInt`](@ref): $(PDATA.iuse_viscous_strain_soft.description)
- `vsoftfac::`[`ParameterFloat`](@ref): $(PDATA.vsoftfac.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of parameter objects

# Nested Dot Access
- `iuse_viscous_strain_soft = $(ROOT_NAME).$(GRP_NAME).iuse_viscous_strain_soft.value`
- `vsoftfac = $(ROOT_NAME).$(GRP_NAME).vsoftfac.value`

# Constructor
    Softening()

Create a new Softening parameter group with default values.

# Returns
- `Softening`: New Softening parameter group with initialized values

"""
mutable struct Softening <: AbstractParameterGroup
    iuse_viscous_strain_soft::ParameterInt
    vsoftfac::ParameterFloat
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function Softening()::Softening
    pdata = get_eb_parameters()
    data = Softening(
        pdata.iuse_viscous_strain_soft,
        pdata.vsoftfac,
        Union{ParameterFloat, ParameterInt}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
