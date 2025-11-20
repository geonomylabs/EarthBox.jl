module JLDData

using JLD2
import EarthBox.JLDTools: intstr
import ..DataNames: MarkerDataNames
import ..DataNames: get_list
import ...PlotParametersManager: PlotParameters
import ...PlotParametersManager.PlotConversionManager: convert_time_units
import ...PlotParametersManager.PlotConversionManager: convert_length_units
import ...PlotParametersManager.PlotConversionManager: convert_length_array_units

function get_jld_field_filename(
    parameters::PlotParameters
)::String
    mainpath = parameters.paths.mainpath
    ioutput = parameters.time.ioutput
    jld_filename = joinpath(mainpath, "fields_$(intstr(ioutput)).jld")
    return jld_filename
end

function get_jld_topo_filename(
    parameters::PlotParameters
)::String
    mainpath = parameters.paths.mainpath
    ioutput = parameters.time.ioutput
    jld_filename = joinpath(mainpath, "topo_$(intstr(ioutput)).jld")
    return jld_filename
end

""" Get marker data from jld file.

Inputs
------
- parameters::PlotParameters: Plot parameters
- marker_data_names::MarkerDataNames: Marker data names

Returns
-------
- model_time::Float64: Model time
- y_sealevel::Float64: Sea level in plot units
- base_level_shift::Float64: Base level shift in plot units
- marker_data_dict::Dict{String, Vector{Float64}}: Marker data dictionary
"""
function get_jld_marker_data(
    parameters::PlotParameters,
    marker_data_names::MarkerDataNames
)::Tuple{Float64, Float64, Float64, Dict{String, Vector{Float64}}}
    
    jld_filename = get_jld_marker_filename(parameters)
    
    model_time = nothing
    marker_data_dict = Dict{String, Vector{Float64}}()
    y_sealevel = 0.0
    base_level_shift = 0.0

    if !isfile(jld_filename)
        error("JLD file not found: $jld_filename")
    end

    jldopen(jld_filename, "r") do jldfile
        
        # Get sea level data
        y_sealevel = get(jldfile, "y_sealevel", nothing)
        if y_sealevel === nothing
            println("!!! WARNING !!! y_sealevel not found in jld file")
            y_sealevel = 0.0
        end
        y_sealevel = convert_length_units(parameters.conversion, "m", y_sealevel)
        
        # Get base level shift data
        base_level_shift = get(jldfile, "base_level_shift", nothing)
        if base_level_shift === nothing
            println("!!! WARNING !!! base_level_shift not found in jld file")
            base_level_shift = 0.0
        end
        base_level_shift = convert_length_units(
            parameters.conversion, "m", base_level_shift)
        
        # Process marker data
        ifind_time = 0
        for dataname in get_list(marker_data_names)
            if haskey(jldfile, dataname)
                model_time = jldfile["time"]
                time_units = jldfile["time_units"]
                dset = jldfile[dataname]
                if ifind_time == 0
                    (
                        model_time, time_units
                    ) = convert_time_units(
                        parameters.conversion, model_time, time_units)
                    ifind_time = 1
                end
                array = copy(dset["array"])
                if dataname in ["marker_x", "marker_y"]
                    length_units = dset["units"]
                    array = convert_length_array_units(
                        parameters.conversion, length_units, array)
                end
            else
                array = zeros(Float64, 1)
            end
            marker_data_dict[dataname] = array
        end
    end
    
    if model_time === nothing
        error("Model time is nothing for jld markers file: $jld_filename")
    end
    
    if isempty(marker_data_dict)
        error("No marker arrays were found in jld marker file: $jld_filename")
    end
    
    return model_time, y_sealevel, base_level_shift, marker_data_dict
end

function get_jld_marker_filename(
    parameters::PlotParameters
)::String
    mainpath = parameters.paths.mainpath
    ioutput = parameters.time.ioutput
    return make_jld_marker_filename(mainpath, ioutput)
end

function make_jld_marker_filename(
    output_dir_path::String,
    ioutput::Int
)::String
    jld_filename = joinpath(output_dir_path, "markers_$(intstr(ioutput)).jld")
    return jld_filename
end

end # module
