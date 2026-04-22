module TimeStep

import CodecZlib: ZlibCompressor
import JLD2
import EarthBox.EarthBoxDtypes: AbstractXDMFTimeStep
import ...OutputDTypes: Markers2djld
import ...OutputDTypes: ScalarFieldMeta
import ...XdmfParts: get_xdmf_start_of_timestep_grid
import ...XdmfParts: get_xdmf_end_of_timestep_grid
import ...XdmfParts: get_xdmf_topology_polyvertex
import ...XdmfParts: get_xdmf_geometry_2d_xy_markers
import ...XdmfParts: get_xdmf_scalar_attribute_on_nodes_for_markers
import ...XdmfUtils: getoutform

struct MarkersXdmfTimeStep <: AbstractXDMFTimeStep
    markers2djld::Markers2djld
    model_time::Float64
    time_units::String
    noutput::Int
    y_sealevel::Float64
    base_level_shift::Float64
    scalar_metas::Vector{ScalarFieldMeta}
    enabled_marker_objs::Vector
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
    
    for scalar_meta in markers_xdmf.scalar_metas
        string *= get_xdmf_scalar_attribute_on_nodes_for_markers(
            markers_xdmf.markers2djld.nmarkers,
            markers_xdmf.markers2djld.jld_markerfile,
            scalar_meta
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
    # Python earthbox markers.py: h5py gzip compression_opts=1 — zlib deflate level 1.
    # Import ZlibCompressor here so CodecZlib is loaded with EarthBox (JLD2's dynamic import
    # can fail when Julia is not started with a project that installs CodecZlib).
    JLD2.jldopen(jld_marker_file_path, "w"; compress=ZlibCompressor(level=1)) do file
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

        # Stream scalar fields one at a time: compute getoutform, write, discard.
        # This keeps at most one nmarkers-sized array in memory at a time instead
        # of accumulating all fields simultaneously before any write occurs.
        for (obj, meta) in zip(markers_xdmf.enabled_marker_objs, markers_xdmf.scalar_metas)
            group = JLD2.Group(file, meta.jld_dataname)
            group["array"] = getoutform(obj)
            group["name"] = meta.name
            group["units"] = meta.units
        end
    end
end

end