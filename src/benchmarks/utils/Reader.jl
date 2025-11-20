module Reader

import EarthBox.PrintFuncs: print_info
import EarthBox.JLDTools: get_jld_marker_id_data
import EarthBox.JLDTools: get_marker_coordinate_arrays
import EarthBox.JLDTools: get_jld_data
import EarthBox.JLDTools: get_jld_topo_data
import EarthBox.JLDTools: intstr
import EarthBox.Arrays.ArrayUtils: check_limits

function read_marker_material_ids(
    input_path::String,
    _base_name::String,
    ioutput::Int
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Int16}, Float64}
    tmyr, marker_matid = get_jld_marker_id_data(ioutput, input_path)
    tmyr, marker_x_m, marker_y_m = get_marker_coordinate_arrays(ioutput, input_path)
    return marker_x_m, marker_y_m, marker_matid, tmyr
end

function read_scalar_grid(
    input_path::String,
    base_name::String,
    ioutput::Int
)::Tuple{Vector{Float64}, Vector{Float64}, Matrix{Float64}, Float64}
    jld_filename = joinpath(input_path, "fields_" * intstr(ioutput) * ".jld")
    print_info("Reading jld file $jld_filename", level=1)
    tmyr, gridy, gridx, vals, _units_dict = get_jld_data(base_name, jld_filename)
    return gridx, gridy, vals, tmyr
end

function read_topography(
    input_path::String,
    _base_name::String,
    ioutput::Int
)::Tuple{Vector{Float64}, Vector{Float64}}
    jld_filename = joinpath(input_path, "topo_" * intstr(ioutput) * ".jld")
    print_info("Reading jld file $jld_filename", level=1)
    topoy, topox, _units = get_jld_topo_data(jld_filename)
    return topox, topoy
end

function read_vector_grid(
    input_path::String,
    _base_name::String,
    ioutput::Int
)::Tuple{Vector{Float64}, Vector{Float64}, Matrix{Float64}, Matrix{Float64}, Float64}
    jld_filename = joinpath(input_path, "vel_cmyr_" * intstr(ioutput) * ".jld")
    print_info("Reading jld file $jld_filename", level=1)
    tmyr, gridy, gridx, vx, _units_dict = get_jld_data("velx_cmyr", jld_filename)
    _tmyr, _gridy, _gridx, vy, _units_dict = get_jld_data("vely_cmyr", jld_filename)
    check_limits("vx", vx)
    check_limits("vy", vy)
    return gridx, gridy, vx, vy, tmyr
end

end  # end of module 