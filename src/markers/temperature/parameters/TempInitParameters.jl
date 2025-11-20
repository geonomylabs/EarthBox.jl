module TempInitParameters

import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: set_parameters!
import EarthBox.ModelDataContainer: set_parameter!
import EarthBox.ParameterGroupTools: set_group_parameters!
import EarthBox.DictUtils: check_dictionary_names
import ..Options: option_names

function load_parameters!(
    model::ModelData,
    option_symbol::Union{Symbol, Nothing};
    kwargs...
)::Nothing
    # Set dictionary parameters
    parameters = get(kwargs, :parameters, nothing)
    check_keys_in_parameters_dict(parameters, option_symbol)
    set_parameters!(model, parameters)
    # Set individual parameters
    adiabatic_gradient = get(kwargs, :adiabatic_gradient, nothing)
    temperature_uniform = get(kwargs, :temperature_uniform, nothing)
    set_parameter!(model, "adiabatic_gradient", adiabatic_gradient)
    set_parameter!(model, "temperature_uniform", temperature_uniform)
    return nothing
end

""" Check if the keys in the parameters dictionary are valid for the given option symbol.
"""
function check_keys_in_parameters_dict(
    parameters::Union{Dict{String, <:Union{Float64, Int, String}}, Nothing},
    option_symbol::Union{Symbol, Nothing}
)::Nothing
    if parameters !== nothing && option_symbol !== nothing
        check_parameters(parameters, Val(option_symbol))
    end
    return nothing
end

Base.@kwdef struct FractureZoneKeyNames
    temperature_top::String = "temperature_top"
    temperature_bottom::String = "temperature_bottom"
    age_lithosphere_left::String = "age_lithosphere_left"
    age_lithosphere_right::String = "age_lithosphere_right"
    adiabatic_gradient::String = "adiabatic_gradient"
    thermal_lithosphere_depth_left::String = "thermal_lithosphere_depth_left"
    thermal_lithosphere_depth_right::String = "thermal_lithosphere_depth_right"
    thermal_diffusivity::String = "thermal_diffusivity"
end

function check_parameters(
    parameters::Union{Dict{String, <:Union{Float64, Int, String}}, Nothing},
    ::Val{option_names.FractureZone}
)::Nothing
    if parameters !== nothing
        key_names = FractureZoneKeyNames()
        valid_keys = [
            key_names.temperature_top,
            key_names.temperature_bottom,
            key_names.age_lithosphere_left,
            key_names.age_lithosphere_right,
            key_names.adiabatic_gradient,
            key_names.thermal_lithosphere_depth_left,
            key_names.thermal_lithosphere_depth_right,
            key_names.thermal_diffusivity
        ]
        check_dictionary_names(valid_keys, parameters)
    end
    return nothing
end

Base.@kwdef struct TemperatureWaveKeyNames
    temperature_top::String = "temperature_top"
    temperature_bottom::String = "temperature_bottom"
    temperature_of_wave::String = "temperature_of_wave"
end

function check_parameters(
    parameters::Union{Dict{String, <:Union{Float64, Int, String}}, Nothing},
    ::Val{option_names.TemperatureWave}
)::Nothing
    if parameters !== nothing
        key_names = TemperatureWaveKeyNames()
        valid_keys = [
            key_names.temperature_top,
            key_names.temperature_bottom,
            key_names.temperature_of_wave
        ]
        check_dictionary_names(valid_keys, parameters)
    end
    return nothing
end

Base.@kwdef struct BoxConvectionKeyNames
    temperature_top::String = "temperature_top"
    temperature_bottom::String = "temperature_bottom"
end

function check_parameters(
    parameters::Union{Dict{String, <:Union{Float64, Int, String}}, Nothing},
    ::Val{option_names.BoxConvection}
)::Nothing
    if parameters !== nothing
        key_names = BoxConvectionKeyNames()
        valid_keys = [
            key_names.temperature_top,
            key_names.temperature_bottom
        ]
        check_dictionary_names(valid_keys, parameters)
    end
    return nothing
end

Base.@kwdef struct LinearKeyNames
    temperature_top_kelvins::String = "temperature_top_kelvins"
    temperature_bottom_kelvins::String = "temperature_bottom_kelvins"
end

function check_parameters(
    parameters::Union{Dict{String, <:Union{Float64, Int, String}}, Nothing},
    ::Val{option_names.Linear}
)::Nothing
    if parameters !== nothing
        key_names = LinearKeyNames()
        valid_keys = [
            key_names.temperature_top_kelvins,
            key_names.temperature_bottom_kelvins
        ]
        check_dictionary_names(valid_keys, parameters)
    end
    return nothing
end

Base.@kwdef struct FourLinearSegmentsKeyNames
    amplitude_perturbation::String = "amplitude_perturbation"
    width_perturbation::String = "width_perturbation"
    temperature_surface::String = "temperature_surface"
    temperature_moho::String = "temperature_moho"
    temperature_base_lith::String = "temperature_base_lith"
    temperature_bottom::String = "temperature_bottom"
end

function check_parameters(
    parameters::Union{Dict{String, <:Union{Float64, Int, String}}, Nothing},
    ::Val{option_names.FourLinearSegments}
)::Nothing
    if parameters !== nothing
        key_names = FourLinearSegmentsKeyNames()
        valid_keys = [
            key_names.amplitude_perturbation,
            key_names.width_perturbation,
            key_names.temperature_surface,
            key_names.temperature_moho,
            key_names.temperature_base_lith,
            key_names.temperature_bottom
        ]
        check_dictionary_names(valid_keys, parameters)
    end
    return nothing
end

Base.@kwdef struct AnalyticalThreeLayerKeyNames
    temperature_top::String = "temperature_top"
    temperature_bottom::String = "temperature_bottom"
    thick_thermal_lithosphere::String = "thick_thermal_lithosphere"
    amplitude_perturbation::String = "amplitude_perturbation"
    width_perturbation::String = "width_perturbation"
    conductivity_upper_crust::String = "conductivity_upper_crust"
    conductivity_lower_crust::String = "conductivity_lower_crust"
    conductivity_mantle::String = "conductivity_mantle"
    heat_production_upper_crust::String = "heat_production_upper_crust"
    heat_production_lower_crust::String = "heat_production_lower_crust"
    heat_production_mantle::String = "heat_production_mantle"
    temperature_base_lith::String = "temperature_base_lith"
    adiabatic_gradient::String = "adiabatic_gradient"
end

function check_parameters(
    parameters::Union{Dict{String, <:Union{Float64, Int, String}}, Nothing},
    ::Val{option_names.AnalyticalThreeLayer}
)::Nothing
    if parameters !== nothing
        key_names = AnalyticalThreeLayerKeyNames()
        valid_keys = [
            key_names.temperature_top,
            key_names.temperature_bottom,
            key_names.thick_thermal_lithosphere,
            key_names.amplitude_perturbation,
            key_names.width_perturbation,
            key_names.conductivity_upper_crust,
            key_names.conductivity_lower_crust,
            key_names.conductivity_mantle,
            key_names.heat_production_upper_crust,
            key_names.heat_production_lower_crust,
            key_names.heat_production_mantle,
            key_names.temperature_base_lith,
            key_names.adiabatic_gradient
        ]
        check_dictionary_names(valid_keys, parameters)
    end
    return nothing
end

function check_for_missing_key_names(
    parameters::Dict{String, <:Union{Float64, Int}},
    valid_keys::Vector{String}
)::Nothing
    for key in valid_keys
        if !haskey(parameters, key)
            error("Missing key name: $key. Add this key to the input parameters dictionary.")
        end
    end
    return nothing
end


end # module 