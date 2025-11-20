module PlotDict

import EarthBox.PrintFuncs: print_list
import ..ReadInput: get_plot_dict
import ..PlotDtypes: PlotDictType, PlotDictTemplateType

""" Struct for scalar plot parameter group names.

The names must match the keys in the plot_dict_template.yml file. The file
is located in the same directory as this file:

    src/earthbox/plot_tools/model_plots/utils/plot_dict_template.py
"""
Base.@kwdef struct ScalarPlotParameterGroupNames
    thermal_conductivity_plot::String = "thermal_conductivity_plot"
    density_plot::String = "density_plot"
    temperature_plot::String = "temperature_plot"
    viscosity_plot::String = "viscosity_plot"
    strain_rate_plot::String = "strain_rate_plot"
    pressure_plot::String = "pressure_plot"
    normal_stress_plot::String = "normal_stress_plot"
    shear_stress_plot::String = "shear_stress_plot"
    shear_plastic_failure_plot::String = "shear_plastic_failure_plot"
    normal_plastic_failure_plot::String = "normal_plastic_failure_plot"
    velocity_plot::String = "velocity_plot"
end

""" See plot_dict_template.yml for details.
"""
function get_valid_parameter_group_keys()::Vector{String}
    template = get_plot_dict_template()
    return collect(keys(template))
end

function get_plot_dict_template()::PlotDictType
    current_module_dir_path = dirname(@__FILE__)
    template_path = joinpath(current_module_dir_path, "plot_dict_template.yml")
    return get_plot_dict(template_path)
end

function copy_template_to_plot_dict(
    plot_dict_template::PlotDictTemplateType
)::PlotDictType
    plot_dict = Dict{String,Any}()
    for (group_key, group_dict) in plot_dict_template
        plot_dict[group_key] = Dict{String,Any}()
        for (param_name, param_list) in group_dict
            plot_dict[group_key][param_name] = param_list[1]
        end
    end
    return plot_dict
end

""" Update general parameters in plot_dict with input values.

# Arguments
- `group_key::String`: Key for the parameter group to update
- `plot_dict::PlotDictType`: Dictionary containing plot parameters
- `kwargs...`: Keyword arguments containing parameter updates
"""
function update_plot_dict!(
    group_key::String,
    plot_dict::PlotDictType;
    kwargs...
)::Nothing
    valid_group_keys = get_valid_parameter_group_keys()
    if group_key in valid_group_keys
        parameter_keys = keys(plot_dict[group_key])
        for (keyword, value) in kwargs
            keyword_string = String(keyword)
            if keyword_string in parameter_keys && !isnothing(value)
                if value isa Symbol
                    value = string(value)
                end
                plot_dict[group_key][keyword_string] = value
            end
            # Commented out to allow for more flexible input
            #else
            #    print_list("Valid Parameters for $group_key", collect(parameter_keys))
            #    error("$keyword_string is not a valid parameter key for $group_key. " *
            #          "Valid keys printed above.")
            #end
        end
    else
        print_list("Valid Groups", collect(valid_group_keys))
        error("$group_key is not a valid parameter type key. " *
              "Valid keys are printed above.")
    end
    return nothing
end

""" Check all keys and types in plot_dict against template.

Types are compared for the first list element, which is the value
of a given parameter.
"""
function check_plot_dict(
    plot_dict::PlotDictType,
    plot_dict_template::PlotDictType
)::Nothing
    for (group_key, group_dict) in plot_dict
        if !(group_key in keys(plot_dict_template))
            error("$group_key is not valid. Check input file plots.yaml")
        end
        for (param_name, param_value) in group_dict
            if !(param_name in keys(plot_dict_template[group_key]))
                error("$param_name is not a valid parameter name. Check inputs.")
            end
            param_value_type = typeof(param_value)

            template_value = plot_dict_template[group_key][param_name]
            if isa(template_value, Char)
                template_value = string(template_value)
            end
            template_value_type = typeof(template_value)

            if param_value_type != template_value_type
                error("Parameter $param_name with value $param_value " *
                      "has type $param_value_type. " *
                      "However, the type should be $template_value_type. " *
                      "Check model plots inputs")
            end
        end
    end
    return nothing
end

end # module 