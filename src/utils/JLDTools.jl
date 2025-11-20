module JLDTools

using JLD2
using Printf
import EarthBox.Arrays.ArrayUtils: check_limits

""" Convert integer to string with leading zeros.
"""
function intstr(number::Int)::String
    return @sprintf("%05d", number)
end

""" Get topography data from jld file.

Returns
-------
- topoy_m::Vector{Float64}: Y topography coordinates in meters
- topox_m::Vector{Float64}: X topography coordinates in meters
- units::String: Units of the coordinates
"""
function get_jld_topo_data(jld_filename::String)::Tuple{Vector{Float64}, Vector{Float64}, String}
    jldopen(jld_filename, "r") do file
        topo_xy_km = file["topo_xy_km"]["array"]
        units = file["topo_xy_km"]["units"]
        topoy_m, topox_m = get_topo_component_arrays(topo_xy_km)
        units = "m"
        return topoy_m, topox_m, units
    end
end

""" Flip y-axis since array is coming from Paraview format.
"""
function get_topo_component_arrays(
    topo_xy_km::Matrix{Float64}
)::Tuple{Vector{Float64}, Vector{Float64}}
    ntopo = size(topo_xy_km, 1)
    topox_m = zeros(Float64, ntopo)
    topoy_m = zeros(Float64, ntopo)
    for i in 1:ntopo
        topox_m[i] = topo_xy_km[i, 1] * 1000.0
        topoy_m[i] = -topo_xy_km[i, 2] * 1000.0
    end
    return topoy_m, topox_m
end

""" Get grid data from jld file.

Returns
-------
- model_time::Float64: Model time
- gridy_km::Vector{Float64}: Y grid coordinates in km
- gridx_km::Vector{Float64}: X grid coordinates in km
- scalar_array::Matrix{Float64}: Scalar array
- units_dict::Dict{String, String}:
    - "length_units": Length units
    - "scalar_units": Scalar units
    - "time_units": Time units
"""
function get_jld_data(
    dataname::String, 
    jld_filename::String;
    print_all_data::Bool=true
)::Tuple{Float64, Vector{Float64}, Vector{Float64}, Matrix{Float64}, Dict{String, String}}
    jldopen(jld_filename, "r") do file
        group = file[dataname]
        model_time = group["model_time"]
        time_units = group["time_units"]
        scalar_array = group["array"]
        grid_type = group["grid_type"]
        gridy_km, gridx_km, length_units = get_jld_grid_data(file, grid_type)
        units_dict = Dict(
            "length_units" => length_units,
            "scalar_units" => group["units"],
            "time_units" => time_units
        )
        return model_time, gridy_km, gridx_km, scalar_array, units_dict
    end
end

""" Get grid data from jld file based on grid type.
"""
function get_jld_grid_data(
    file::JLD2.JLDFile, 
    grid_type::String
)::Tuple{Vector{Float64}, Vector{Float64}, String}
    if grid_type == "basic"
        gridy_km = -file["gridy_km"]["array"]
        units = file["gridy_km"]["units"]
        gridx_km = file["gridx_km"]["array"]
    elseif grid_type == "pressure"
        gridy_km = -file["gridy_pr_km"]["array"]
        units = file["gridy_pr_km"]["units"]
        gridx_km = file["gridx_pr_km"]["array"]
    else
        error("Unknown grid type: $grid_type")
    end
    return gridy_km, gridx_km, units
end

""" Get marker coordinates from output file.

Inputs
------
- itime_step::Int: Time step integer id
- dir_path::String: Path to the directory where the jld file is located

Returns
-------
- model_time::Float64: Model time (Myr)
- marker_x_m::Vector{Float64}: Marker x coordinates in meters
- marker_y_m::Vector{Float64}: Marker y coordinates in meters
"""
function get_marker_coordinate_arrays(
    itime_step::Int, 
    dir_path::String
)::Tuple{Float64, Vector{Float64}, Vector{Float64}}
    jld_filename = get_jld_marker_filename(itime_step, dir_path)
    jldopen(jld_filename, "r") do file
        marker_x_m = file["marker_x"]["array"]
        marker_y_m = file["marker_y"]["array"]
        model_time = file["time"]
        return model_time, marker_x_m, marker_y_m
    end
end

""" Get marker material id's from jld file.

Inputs
------
- itime_step::Int: Time step integer id
- dir_path::String: Path to the directory where the jld file is located

Returns
-------
- model_time::Float64: Model time (Myr)
- marker_matid::Vector{Int16}: Marker material id's
"""
function get_jld_marker_id_data(
    itime_step::Int, 
    dir_path::String
)::Tuple{Float64, Vector{Int16}}
    jld_filename = get_jld_marker_filename(itime_step, dir_path)
    jldopen(jld_filename, "r") do file
        marker_matid = convert(Vector{Int16}, file["marker_matid"]["array"])
        model_time = file["time"]
        return model_time, marker_matid
    end
end

""" Get marker data from jld file.

Inputs
------
- itime_step::Int: Time step integer id
- dir_path::String: Path to the directory where the jld file is located
- dataname::String: Name of the data to retrieve

Returns
-------
- model_time::Float64: Model time (Myr)
- marker_data::Vector{Float64}: Marker data from jld file
"""
function get_jld_marker_data(
    itime_step::Int, 
    dir_path::String, 
    dataname::String
)::Tuple{Float64, Vector{Float64}}
    jld_filename = get_jld_marker_filename(itime_step, dir_path)
    jldopen(jld_filename, "r") do file
        marker_data = file[dataname]["array"]
        model_time = file["time"]
        return model_time, marker_data
    end
end

""" Get the jld marker filename.

Inputs
------
- itime_step::Int: Time step integer id
- dir_path::String: Path to the directory where the jld file is located
"""
function get_jld_marker_filename(itime_step::Int, dir_path::String)::String
    return joinpath(dir_path, "markers_" * intstr(itime_step) * ".jld")
end

""" Get marker data from jld file.

Inputs
------
- ioutput::Int: Output integer id
- output_dir_path::String: Path to the directory where the jld file is located

Returns
-------
- y_sealevel::Float64: Y-coordinate of sea level (meters)
"""
function get_y_sealevel_from_jld_marker_data(
    ioutput::Int, 
    output_dir_path::String
)::Float64
    jld_filename = get_jld_marker_filename(ioutput, output_dir_path)
    jldopen(jld_filename, "r") do file
        return file["y_sealevel"]
    end
end

""" Get marker data from jld file.

Inputs
------
- ioutput::Int: Output integer id
- output_dir_path::String: Path to the directory where the jld file is located

Returns
-------
- base_level_shift::Float64: Base level shift (meters) (Positive = downward shift)
"""
function get_base_level_shift_from_jld_marker_data(
    ioutput::Int, 
    output_dir_path::String
)::Float64
    jld_filename = get_jld_marker_filename(ioutput, output_dir_path)
    jldopen(jld_filename, "r") do file
        base_level_shift = file["base_level_shift"]
    end
    if ioutput == 1
        base_level_shift = 0.0
    end
    return base_level_shift
end 

function print_all_data_names(
    jld_filename::String,
    file::JLD2.JLDFile,
    print_all_data::Bool=true
)::Nothing
    if print_all_data
        println("Available data in $jld_filename:")
        for key in keys(file)
            println("  - $key")
            if typeof(file[key]) <: JLD2.Group
                println("    Group contents:")
                for subkey in keys(file[key])
                    println("      - $subkey")
                end
            end
        end
    end
    return nothing
end


end # module
