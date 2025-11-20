module TimeStep

import EarthBox.EarthBoxDtypes: AbstractXDMFTimeStep
import JLD2
import ...OutputDTypes: Topo2djld
import ...XdmfParts: get_xdmf_start_of_timestep_grid
import ...XdmfParts: get_xdmf_end_of_timestep_grid
import ...XdmfParts: get_xdmf_topology_polyvertex
import ...XdmfParts: get_xdmf_geometry_2d_xy_markers
import ...XdmfTimeStepTools: set_time, define_time_attr

struct TopographyXdmfTimeStep <: AbstractXDMFTimeStep
    topo2djld::Topo2djld
    model_time::Float64
    time_units::String
    noutput::Int
end

function get_xdmf_string_for_timestep(step::TopographyXdmfTimeStep)
    string = get_xdmf_start_of_timestep_grid(
        "topo", step.model_time, step.time_units)
    string *= get_xdmf_topology_polyvertex(
        step.topo2djld.ntopo, step.topo2djld.jld_markerfile,
        step.topo2djld.jld_dataname_ids)
    string *= get_xdmf_geometry_2d_xy_markers(
        step.topo2djld.ntopo, step.topo2djld.jld_markerfile,
        step.topo2djld.jld_dataname_xy)
    string *= get_xdmf_end_of_timestep_grid()
    return string
end

function make_jld2_file(data::TopographyXdmfTimeStep, output_dir::String)
    jld_marker_filename = data.topo2djld.jld_markerfile
    jld_dataname_xy = data.topo2djld.jld_dataname_xy
    jld_dataname_ids = data.topo2djld.jld_dataname_ids
    topo_id_array = data.topo2djld.topo_id_array
    topo_xy_km_array = data.topo2djld.topo_xy_km_array

    jld_marker_file_path = joinpath(output_dir, jld_marker_filename)
    
    JLD2.jldopen(jld_marker_file_path, "w") do file
        set_time(data, file)
        
        # Save marker IDs
        group = JLD2.Group(file, jld_dataname_ids)
        group["array"] = topo_id_array
        group["name"] = "marker_ids"
        group["units"] = "None"
        define_time_attr(data, group)
        
        # Save marker XY coordinates
        group = JLD2.Group(file, jld_dataname_xy)
        group["array"] = topo_xy_km_array
        group["name"] = "marker_xy"
        group["units"] = "km"
        define_time_attr(data, group)
    end
end

end # module
