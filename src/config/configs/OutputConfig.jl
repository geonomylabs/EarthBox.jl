module OutputConfig

import EarthBox.PrintFuncs: print_info
import YAML

Base.@kwdef struct MarkerOutputKeys
    marker_x::String = "marker_x"
    marker_y::String = "marker_y"
    marker_matid::String = "marker_matid"
    marker_GII::String = "marker_GII"
    marker_strain_plastic::String = "marker_strain_plastic"
    marker_strain_rate_plastic::String = "marker_strain_rate_plastic"
    marker_meltfrac::String = "marker_meltfrac"
    marker_extractable_meltfrac::String = "marker_extractable_meltfrac"
    marker_extracted_meltfrac::String = "marker_extracted_meltfrac"
    marker_serpentinization::String = "marker_serpentinization"
    marker_pfailure::String = "marker_pfailure"
    marker_TK::String = "marker_TK"
    marker_xn::String = "marker_xn"
    marker_yn::String = "marker_yn"
    marker_sxx::String = "marker_sxx"
    marker_sxy::String = "marker_sxy"
    marker_eta::String = "marker_eta"
    marker_exx::String = "marker_exx"
    marker_exy::String = "marker_exy"
    marker_pr::String = "marker_pr"
    marker_sr_ratio::String = "marker_sr_ratio"
    marker_age::String = "marker_age"
    marker_fric_ini::String = "marker_fric_ini"
    marker_fric::String = "marker_fric"
    marker_cohesion::String = "marker_cohesion"
    marker_preexp::String = "marker_preexp"
    marker_rho::String = "marker_rho"
end

"""
    OutputConfigState

Mutable struct to store output configuration parameters.

# Fields
- `general::Dict{String, Bool}`: General output configuration
- `marker_output::Dict{String, Bool}`: Marker output configuration
- `marker_output_keys::MarkerOutputKeys`: Marker output keys
"""
mutable struct OutputConfigState
    general::Dict{String, Bool}
    marker_output::Dict{String, Bool}
    marker_output_keys::MarkerOutputKeys
end

function OutputConfigState()
    marker_keys = MarkerOutputKeys()
    general = Dict{String, Bool}("make_files" => true)
    marker_output = Dict{String, Bool}(
        marker_keys.marker_x                    => true,
        marker_keys.marker_y                    => true,
        marker_keys.marker_matid                => true,
        marker_keys.marker_GII                  => true,
        marker_keys.marker_strain_plastic       => true,
        marker_keys.marker_strain_rate_plastic  => true,
        marker_keys.marker_meltfrac             => true,
        marker_keys.marker_extractable_meltfrac => false,
        marker_keys.marker_extracted_meltfrac   => false,
        marker_keys.marker_serpentinization     => false,
        marker_keys.marker_pfailure             => true,
        marker_keys.marker_TK                   => false,
        marker_keys.marker_xn                   => false,
        marker_keys.marker_yn                   => false,
        marker_keys.marker_sxx                  => false,
        marker_keys.marker_sxy                  => false,
        marker_keys.marker_eta                  => false,
        marker_keys.marker_exx                  => false,
        marker_keys.marker_exy                  => false,
        marker_keys.marker_pr                   => false,
        marker_keys.marker_sr_ratio             => false,
        marker_keys.marker_age                  => true,
        marker_keys.marker_fric_ini             => false,
        marker_keys.marker_fric                 => false,
        marker_keys.marker_cohesion             => false,
        marker_keys.marker_preexp               => false,
        marker_keys.marker_rho                  => false
    )
    return OutputConfigState(general, marker_output, marker_keys)
end

function read_marker_output_config!(
    config::OutputConfigState, 
    output_options_file::String
)::Nothing
    output_config = YAML.load_file(output_options_file)
    check_and_update_output_config!(config, output_config)
    return nothing
end

function check_and_update_output_config!(
    config::OutputConfigState, 
    output_config::Dict
)::Nothing
    (
        general_keys_master, marker_output_keys_master
    ) = get_marker_output_config_master_keys(config)

    general_keys_form_input = collect(keys(output_config["general"]))
    marker_output_keys_form_input = collect(keys(output_config["marker_output"]))

    for master_key in marker_output_keys_master
        if !(master_key in marker_output_keys_form_input)
            println("    >> Warning found missing output config key: ", master_key)
            output_config["marker_output"][master_key] = true
        end
    end

    for master_key in general_keys_master
        if !(master_key in general_keys_form_input)
            println("    >> Warning found missing output config key: ", master_key)
            output_config["general"][master_key] = true
        end
    end
    return nothing
end

function get_marker_output_config_master_keys(config::OutputConfigState)
    general_keys = collect(keys(config.general))
    marker_output_keys = collect(keys(config.marker_output))
    return general_keys, marker_output_keys
end

function print_output_config(config::OutputConfigState)
    print_info("", level=1)
    print_info("Output Configuration", level=1)
    print_info("=================================", level=1)
    print_info("General Output", level=1)
    print_info("---------------------------------", level=1)
    for (k, v) in config.general
        print_info("$k : $v", level=2)
    end
    print_info("---------------------------------", level=1)
    print_info("Marker Output", level=1)
    print_info("---------------------------------", level=1)
    for (k, v) in config.marker_output
        print_info("$k : $v", level=2)
    end
    print_info("", level=1)
end

end # module OutputConfig 