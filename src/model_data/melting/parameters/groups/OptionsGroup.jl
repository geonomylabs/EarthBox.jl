module OptionsGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterInt
import EarthBox.Parameters: print_parameter
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.melting.parameters"
const GRP_NAME = "options"

const PDATA = get_eb_parameters()

"""
    Options <: AbstractParameterGroup

Parameter group for melting model on/off switches.

# Fields
- `iuse_melting::`[`ParameterInt`](@ref): $(PDATA.iuse_melting.description)
- `iuse_melt_viscosity::`[`ParameterInt`](@ref): $(PDATA.iuse_melt_viscosity.description)
- `iuse_melt_thermal_props::`[`ParameterInt`](@ref): $(PDATA.iuse_melt_thermal_props.description)
- `iuse_extraction::`[`ParameterInt`](@ref): $(PDATA.iuse_extraction.description)
- `iuse_gabbroic_fractionation::`[`ParameterInt`](@ref): $(PDATA.iuse_gabbroic_fractionation.description)
- `iuse_shallow_mantle_injection::`[`ParameterInt`](@ref): $(PDATA.iuse_shallow_mantle_injection.description)
- `iuse_random_injection_subdomain::`[`ParameterInt`](@ref): $(PDATA.iuse_random_injection_subdomain.description)
- `iuse_normal_injection_subdomain::`[`ParameterInt`](@ref): $(PDATA.iuse_normal_injection_subdomain.description)
- `iuse_depletion_density::`[`ParameterInt`](@ref): $(PDATA.iuse_depletion_density.description)
- `iuse_exponential_viscosity_reduction::`[`ParameterInt`](@ref): $(PDATA.iuse_exponential_viscosity_reduction.description)
- `obj_list::Vector{ParameterInt}`: List of parameter objects

# Nested Dot Access
- `iuse_melting = $(ROOT_NAME).$(GRP_NAME).iuse_melting.value`
- `iuse_melt_viscosity = $(ROOT_NAME).$(GRP_NAME).iuse_melt_viscosity.value`
- `iuse_melt_thermal_props = $(ROOT_NAME).$(GRP_NAME).iuse_melt_thermal_props.value`
- `iuse_extraction = $(ROOT_NAME).$(GRP_NAME).iuse_extraction.value`
- `iuse_gabbroic_fractionation = $(ROOT_NAME).$(GRP_NAME).iuse_gabbroic_fractionation.value`
- `iuse_shallow_mantle_injection = $(ROOT_NAME).$(GRP_NAME).iuse_shallow_mantle_injection.value`
- `iuse_random_injection_subdomain = $(ROOT_NAME).$(GRP_NAME).iuse_random_injection_subdomain.value`
- `iuse_normal_injection_subdomain = $(ROOT_NAME).$(GRP_NAME).iuse_normal_injection_subdomain.value`
- `iuse_depletion_density = $(ROOT_NAME).$(GRP_NAME).iuse_depletion_density.value`
- `iuse_exponential_viscosity_reduction = $(ROOT_NAME).$(GRP_NAME).iuse_exponential_viscosity_reduction.value`

# Constructor
    Options()

# Returns
- `Options`: New Options parameter group with initialized values

"""
mutable struct Options <: AbstractParameterGroup
    iuse_melting::ParameterInt
    iuse_melt_viscosity::ParameterInt
    iuse_melt_thermal_props::ParameterInt
    iuse_extraction::ParameterInt
    iuse_gabbroic_fractionation::ParameterInt
    iuse_shallow_mantle_injection::ParameterInt
    iuse_random_injection_subdomain::ParameterInt
    iuse_normal_injection_subdomain::ParameterInt
    iuse_depletion_density::ParameterInt
    iuse_exponential_viscosity_reduction::ParameterInt
    obj_list::Vector{ParameterInt}
end

function Options()::Options
    pdata = get_eb_parameters()
    data = Options(
        pdata.iuse_melting,
        pdata.iuse_melt_viscosity,
        pdata.iuse_melt_thermal_props,
        pdata.iuse_extraction,
        pdata.iuse_gabbroic_fractionation,
        pdata.iuse_shallow_mantle_injection,
        pdata.iuse_random_injection_subdomain,
        pdata.iuse_normal_injection_subdomain,
        pdata.iuse_depletion_density,
        pdata.iuse_exponential_viscosity_reduction,
        ParameterInt[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
