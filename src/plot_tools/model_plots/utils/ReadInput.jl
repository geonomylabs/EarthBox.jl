module ReadInput

using YAML
using ..PlotDtypes: PlotDictType, PlotParametersType

function get_plot_dict(input_file_name::String)::PlotDictType
    data_dict = read_yaml_file(input_file_name)
    plot_dict = PlotDictType()
    for (key, subdict) in data_dict
        current_keys = keys(plot_dict)
        if !(key in current_keys)
            plot_dict[key] = PlotParametersType()
            for (subkey, sublist) in subdict
                parameter_value = sublist[1]
                if parameter_value isa Vector{Float64}
                    parameter_tuple = Tuple(parameter_value)
                    plot_dict[key][subkey] = parameter_tuple
                else
                    plot_dict[key][subkey] = parameter_value
                end
            end
        end
    end
    return plot_dict
end

function read_yaml_file(yaml_file_path::String)::Dict{Any,Any}
    return YAML.load_file(yaml_file_path)
end

end # module 