module HeatOptionsGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt, ParameterStr
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.heat_equation.parameters"
const GRP_NAME = "heat_options"

const PDATA = get_eb_parameters()

"""
    HeatOptions <: AbstractParameterGroup

Parameter group for heat solver on/off switches.

# Fields
- `iuse_heat::`[`ParameterInt`](@ref): $(PDATA.iuse_heat.description)
- `iuse_shear_heating::`[`ParameterInt`](@ref): $(PDATA.iuse_shear_heating.description)
- `iuse_adiabatic_heating::`[`ParameterInt`](@ref): $(PDATA.iuse_adiabatic_heating.description)
- `iuse_sticky_correction::`[`ParameterInt`](@ref): $(PDATA.iuse_sticky_correction.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of parameter objects

# Nested Dot Access
- `iuse_heat = $(ROOT_NAME).$(GRP_NAME).iuse_heat.value`
- `iuse_shear_heating = $(ROOT_NAME).$(GRP_NAME).iuse_shear_heating.value`
- `iuse_adiabatic_heating = $(ROOT_NAME).$(GRP_NAME).iuse_adiabatic_heating.value`
- `iuse_sticky_correction = $(ROOT_NAME).$(GRP_NAME).iuse_sticky_correction.value`

# Constructor
    HeatOptions()

# Returns
- `HeatOptions`: New HeatOptions parameter group with initialized values

"""
mutable struct HeatOptions <: AbstractParameterGroup
    iuse_heat::ParameterInt
    iuse_shear_heating::ParameterInt
    iuse_adiabatic_heating::ParameterInt
    iuse_sticky_correction::ParameterInt
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function HeatOptions()::HeatOptions
    pdata = get_eb_parameters()
    data = HeatOptions(
        pdata.iuse_heat,
        pdata.iuse_shear_heating,
        pdata.iuse_adiabatic_heating,
        pdata.iuse_sticky_correction,
        Union{ParameterFloat, ParameterInt}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
