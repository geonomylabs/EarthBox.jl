module InitialConditionGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt, ParameterStr
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.heat_equation.parameters"
const GRP_NAME = "initial_condition"

const PDATA = get_eb_parameters()

"""
    InitialCondition <: AbstractParameterGroup

Parameter group for initial temperature condition parameters.

# Fields
- `itype_temp::`[`ParameterInt`](@ref): $(PDATA.itype_temp.description)
- `stype_temp::`[`ParameterStr`](@ref): $(PDATA.stype_temp.description)
- `temperature_uniform::`[`ParameterFloat`](@ref): $(PDATA.temperature_uniform.description)
- `temperature_of_wave::`[`ParameterFloat`](@ref): $(PDATA.temperature_of_wave.description)
- `age_lithosphere_left::`[`ParameterFloat`](@ref): $(PDATA.age_lithosphere_left.description)
- `age_lithosphere_right::`[`ParameterFloat`](@ref): $(PDATA.age_lithosphere_right.description)
- `thermal_lithosphere_depth_left::`[`ParameterFloat`](@ref): $(PDATA.thermal_lithosphere_depth_left.description)
- `thermal_lithosphere_depth_right::`[`ParameterFloat`](@ref): $(PDATA.thermal_lithosphere_depth_right.description)
- `thermal_diffusivity::`[`ParameterFloat`](@ref): $(PDATA.thermal_diffusivity.description)
- `adiabatic_gradient::`[`ParameterFloat`](@ref): $(PDATA.adiabatic_gradient.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of parameter objects

# Nested Dot Access
- `itype_temp = $(ROOT_NAME).$(GRP_NAME).itype_temp.value`
- `stype_temp = $(ROOT_NAME).$(GRP_NAME).stype_temp.value`
- `temperature_uniform = $(ROOT_NAME).$(GRP_NAME).temperature_uniform.value`
- `temperature_of_wave = $(ROOT_NAME).$(GRP_NAME).temperature_of_wave.value`
- `age_lithosphere_left = $(ROOT_NAME).$(GRP_NAME).age_lithosphere_left.value`
- `age_lithosphere_right = $(ROOT_NAME).$(GRP_NAME).age_lithosphere_right.value`
- `thermal_lithosphere_depth_left = $(ROOT_NAME).$(GRP_NAME).thermal_lithosphere_depth_left.value`
- `thermal_lithosphere_depth_right = $(ROOT_NAME).$(GRP_NAME).thermal_lithosphere_depth_right.value`
- `thermal_diffusivity = $(ROOT_NAME).$(GRP_NAME).thermal_diffusivity.value`
- `adiabatic_gradient = $(ROOT_NAME).$(GRP_NAME).adiabatic_gradient.value`

# Constructor
    InitialCondition()

# Returns
- `InitialCondition`: New InitialCondition parameter group with initialized values

"""
mutable struct InitialCondition <: AbstractParameterGroup
    itype_temp::ParameterInt
    stype_temp::ParameterStr
    temperature_uniform::ParameterFloat
    temperature_of_wave::ParameterFloat
    age_lithosphere_left::ParameterFloat
    age_lithosphere_right::ParameterFloat
    thermal_lithosphere_depth_left::ParameterFloat
    thermal_lithosphere_depth_right::ParameterFloat
    thermal_diffusivity::ParameterFloat
    adiabatic_gradient::ParameterFloat
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function InitialCondition()::InitialCondition
    pdata = get_eb_parameters()
    data = InitialCondition(
        pdata.itype_temp,
        pdata.stype_temp,
        pdata.temperature_uniform,
        pdata.temperature_of_wave,
        pdata.age_lithosphere_left,
        pdata.age_lithosphere_right,
        pdata.thermal_lithosphere_depth_left,
        pdata.thermal_lithosphere_depth_right,
        pdata.thermal_diffusivity,
        pdata.adiabatic_gradient,
        Union{ParameterFloat, ParameterInt}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
