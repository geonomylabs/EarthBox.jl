module TimeStep

import JLD2
import EarthBox.EarthBoxDtypes: AbstractXDMFTimeStep
import ...OutputDTypes: Markers2djld
import ...OutputDTypes: ScalarField
import ...XdmfParts: get_xdmf_start_of_timestep_grid
import ...XdmfParts: get_xdmf_end_of_timestep_grid
import ...XdmfParts: get_xdmf_topology_polyvertex
import ...XdmfParts: get_xdmf_geometry_2d_xy_markers
import ...XdmfParts: get_xdmf_scalar_attribute_on_nodes_for_markers

struct MarkersXdmfTimeStep <: AbstractXDMFTimeStep
    markers2djld::Markers2djld
    model_time::Float64
    time_units::String
    noutput::Int
    y_sealevel::Float64
    base_level_shift::Float64
    scalars_on_markers::Vector{ScalarField}
end

function get_xdmf_string_for_timestep(markers_xdmf::MarkersXdmfTimeStep)::String
    string = get_xdmf_start_of_timestep_grid(
        "points", markers_xdmf.model_time, markers_xdmf.time_units
    )
    string *= get_xdmf_topology_polyvertex(
        markers_xdmf.markers2djld.nmarkers,
        markers_xdmf.markers2djld.jld_markerfile,
        markers_xdmf.markers2djld.jld_dataname_ids
    )
    string *= get_xdmf_geometry_2d_xy_markers(
        markers_xdmf.markers2djld.nmarkers,
        markers_xdmf.markers2djld.jld_markerfile,
        markers_xdmf.markers2djld.jld_dataname_xy
    )
    
    for scalar_field in markers_xdmf.scalars_on_markers
        string *= get_xdmf_scalar_attribute_on_nodes_for_markers(
            markers_xdmf.markers2djld.nmarkers,
            markers_xdmf.markers2djld.jld_markerfile,
            scalar_field
        )
    end
    string *= get_xdmf_end_of_timestep_grid()
    return string
end

function make_jld2_file(markers_xdmf::MarkersXdmfTimeStep, output_dir::String)
    jld_marker_filename = markers_xdmf.markers2djld.jld_markerfile
    jld_dataname_xy = markers_xdmf.markers2djld.jld_dataname_xy
    jld_dataname_ids = markers_xdmf.markers2djld.jld_dataname_ids
    marker_id_array = markers_xdmf.markers2djld.marker_id_array
    marker_xy_km_array = markers_xdmf.markers2djld.marker_xy_km_array

    jld_marker_file_path = joinpath(output_dir, jld_marker_filename)
    
    JLD2.jldopen(jld_marker_file_path, "w") do file
        # Set time attributes
        file["time"] = markers_xdmf.model_time
        file["time_units"] = markers_xdmf.time_units
        file["y_sealevel"] = markers_xdmf.y_sealevel
        file["base_level_shift"] = markers_xdmf.base_level_shift

        # Save marker IDs
        group = JLD2.Group(file, jld_dataname_ids)
        group["array"] = marker_id_array
        group["name"] = "marker_ids"
        group["units"] = "None"

        # Save marker XY coordinates
        group = JLD2.Group(file, jld_dataname_xy)
        group["array"] = marker_xy_km_array
        group["name"] = "marker_xy"
        group["units"] = "km"

        # Save scalar fields
        for scalar_field in markers_xdmf.scalars_on_markers
            jld_dataname = scalar_field.jld_dataname
            scalar_array = scalar_field.scalar_array
            group = JLD2.Group(file, jld_dataname)
            group["array"] = scalar_array
            group["name"] = scalar_field.name
            group["units"] = scalar_field.units
        end
    end
end

end