module ThermalCondGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt, ParameterStr
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.heat_equation.parameters"
const GRP_NAME = "thermalcond"

const PDATA = get_eb_parameters()

"""
    ThermalCond <: AbstractParameterGroup

Parameter group for thermal conductivity model options.

# Fields
- `itype_conductivity::`[`ParameterInt`](@ref): $(PDATA.itype_conductivity.description)
- `stype_conductivity::`[`ParameterStr`](@ref): $(PDATA.stype_conductivity.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of parameter objects

# Nested Dot Access
- `itype_conductivity = $(ROOT_NAME).$(GRP_NAME).itype_conductivity.value`
- `stype_conductivity = $(ROOT_NAME).$(GRP_NAME).stype_conductivity.value`

# Constructor
    ThermalCond()

# Returns
- `ThermalCond`: New ThermalCond parameter group with initialized values

"""
mutable struct ThermalCond <: AbstractParameterGroup
    itype_conductivity::ParameterInt
    stype_conductivity::ParameterStr
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function ThermalCond()::ThermalCond
    pdata = get_eb_parameters()
    data = ThermalCond(
        pdata.itype_conductivity,
        pdata.stype_conductivity,
        Union{ParameterFloat, ParameterInt}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
