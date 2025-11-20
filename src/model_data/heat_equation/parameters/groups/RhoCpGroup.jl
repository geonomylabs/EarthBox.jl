module RhoCpGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt, ParameterStr
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.heat_equation.parameters"
const GRP_NAME = "rhocp"

const PDATA = get_eb_parameters()

"""
    RhoCp <: AbstractParameterGroup

Parameter group for density times heat capacity model options.

# Fields
- `itype_rhocp::`[`ParameterInt`](@ref): $(PDATA.itype_rhocp.description)
- `stype_rhocp::`[`ParameterStr`](@ref): $(PDATA.stype_rhocp.description)
- `maximum_heat_capacity::`[`ParameterFloat`](@ref): $(PDATA.maximum_heat_capacity.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of parameter objects

# Nested Dot Access
- `itype_rhocp = $(ROOT_NAME).$(GRP_NAME).itype_rhocp.value`
- `stype_rhocp = $(ROOT_NAME).$(GRP_NAME).stype_rhocp.value`
- `maximum_heat_capacity = $(ROOT_NAME).$(GRP_NAME).maximum_heat_capacity.value`

# Constructor
    RhoCp()

Create a new RhoCp parameter group with default values.

# Returns
- `RhoCp`: New RhoCp parameter group with initialized values

"""
mutable struct RhoCp <: AbstractParameterGroup
    itype_rhocp::ParameterInt
    stype_rhocp::ParameterStr
    maximum_heat_capacity::ParameterFloat
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function RhoCp()::RhoCp
    pdata = get_eb_parameters()
    data = RhoCp(
        pdata.itype_rhocp,
        pdata.stype_rhocp,
        pdata.maximum_heat_capacity,
        Union{ParameterFloat, ParameterInt}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
