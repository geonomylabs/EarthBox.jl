module SteadyStateGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt, ParameterStr
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.heat_equation.parameters"
const GRP_NAME = "steady_state"

const PDATA = get_eb_parameters()

"""
    SteadyState <: AbstractParameterGroup

Parameter group for steady state thermal structure parameters.

# Fields
- `thick_thermal_lithosphere::`[`ParameterFloat`](@ref): $(PDATA.thick_thermal_lithosphere.description)
- `conductivity_upper_crust::`[`ParameterFloat`](@ref): $(PDATA.conductivity_upper_crust.description)
- `conductivity_lower_crust::`[`ParameterFloat`](@ref): $(PDATA.conductivity_lower_crust.description)
- `conductivity_mantle::`[`ParameterFloat`](@ref): $(PDATA.conductivity_mantle.description)
- `heat_production_upper_crust::`[`ParameterFloat`](@ref): $(PDATA.heat_production_upper_crust.description)
- `heat_production_lower_crust::`[`ParameterFloat`](@ref): $(PDATA.heat_production_lower_crust.description)
- `heat_production_mantle::`[`ParameterFloat`](@ref): $(PDATA.heat_production_mantle.description)
- `amplitude_perturbation::`[`ParameterFloat`](@ref): $(PDATA.amplitude_perturbation.description)
- `width_perturbation::`[`ParameterFloat`](@ref): $(PDATA.width_perturbation.description)
- `temperature_surface::`[`ParameterFloat`](@ref): $(PDATA.temperature_surface.description)
- `temperature_moho::`[`ParameterFloat`](@ref): $(PDATA.temperature_moho.description)
- `temperature_base_lith::`[`ParameterFloat`](@ref): $(PDATA.temperature_base_lith.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of parameter objects

# Nested Dot Access
- `thick_thermal_lithosphere = $(ROOT_NAME).$(GRP_NAME).thick_thermal_lithosphere.value`
- `conductivity_upper_crust = $(ROOT_NAME).$(GRP_NAME).conductivity_upper_crust.value`
- `conductivity_lower_crust = $(ROOT_NAME).$(GRP_NAME).conductivity_lower_crust.value`
- `conductivity_mantle = $(ROOT_NAME).$(GRP_NAME).conductivity_mantle.value`
- `heat_production_upper_crust = $(ROOT_NAME).$(GRP_NAME).heat_production_upper_crust.value`
- `heat_production_lower_crust = $(ROOT_NAME).$(GRP_NAME).heat_production_lower_crust.value`
- `heat_production_mantle = $(ROOT_NAME).$(GRP_NAME).heat_production_mantle.value`
- `amplitude_perturbation = $(ROOT_NAME).$(GRP_NAME).amplitude_perturbation.value`
- `width_perturbation = $(ROOT_NAME).$(GRP_NAME).width_perturbation.value`
- `temperature_surface = $(ROOT_NAME).$(GRP_NAME).temperature_surface.value`
- `temperature_moho = $(ROOT_NAME).$(GRP_NAME).temperature_moho.value`
- `temperature_base_lith = $(ROOT_NAME).$(GRP_NAME).temperature_base_lith.value`

# Constructor
    SteadyState()

# Returns
- `SteadyState`: New SteadyState parameter group with initialized values

"""
mutable struct SteadyState <: AbstractParameterGroup
    thick_thermal_lithosphere::ParameterFloat
    conductivity_upper_crust::ParameterFloat
    conductivity_lower_crust::ParameterFloat
    conductivity_mantle::ParameterFloat
    heat_production_upper_crust::ParameterFloat
    heat_production_lower_crust::ParameterFloat
    heat_production_mantle::ParameterFloat
    amplitude_perturbation::ParameterFloat
    width_perturbation::ParameterFloat
    temperature_surface::ParameterFloat
    temperature_moho::ParameterFloat
    temperature_base_lith::ParameterFloat
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function SteadyState()::SteadyState
    pdata = get_eb_parameters()
    data = SteadyState(
        pdata.thick_thermal_lithosphere,
        pdata.conductivity_upper_crust,
        pdata.conductivity_lower_crust,
        pdata.conductivity_mantle,
        pdata.heat_production_upper_crust,
        pdata.heat_production_lower_crust,
        pdata.heat_production_mantle,
        pdata.amplitude_perturbation,
        pdata.width_perturbation,
        pdata.temperature_surface,
        pdata.temperature_moho,
        pdata.temperature_base_lith,
        Union{ParameterFloat, ParameterInt}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
