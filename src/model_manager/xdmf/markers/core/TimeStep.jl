module TimeStep

import CodecZlib: ZlibCompressor
import JLD2
import EarthBox.EarthBoxDtypes: AbstractXDMFTimeStep
import ...OutputDTypes: Markers2djld
import ...OutputDTypes: ScalarFieldMeta
import ...XdmfParts: get_xdmf_start_of_timestep_grid
import ...XdmfParts: get_xdmf_end_of_timestep_grid
import ...XdmfParts: get_xdmf_topology_polyvertex
import ...XdmfParts: get_xdmf_geometry_2d_vxvy_markers
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
    marker_x_obj
    marker_y_obj
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
    string *= get_xdmf_geometry_2d_vxvy_markers(
        markers_xdmf.markers2djld.nmarkers,
        markers_xdmf.markers2djld.jld_markerfile,
        markers_xdmf.markers2djld.jld_dataname_x,
        markers_xdmf.markers2djld.jld_dataname_y
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
    jld_dataname_x = markers_xdmf.markers2djld.jld_dataname_x
    jld_dataname_y = markers_xdmf.markers2djld.jld_dataname_y
    jld_dataname_ids = markers_xdmf.markers2djld.jld_dataname_ids
    nmarkers = markers_xdmf.markers2djld.nmarkers

    jld_marker_file_path = joinpath(output_dir, jld_marker_filename)
    println(">> RAM make_jld2_file start: $(round(Sys.maxrss()/1024/1024, digits=2)) MB")
    # Python earthbox markers.py: h5py gzip compression_opts=1 — zlib deflate level 1.
    # Import ZlibCompressor here so CodecZlib is loaded with EarthBox (JLD2's dynamic import
    # can fail when Julia is not started with a project that installs CodecZlib).
    JLD2.jldopen(jld_marker_file_path, "w"; compress=ZlibCompressor(level=1)) do file
        file["time"] = markers_xdmf.model_time
        file["time_units"] = markers_xdmf.time_units
        file["y_sealevel"] = markers_xdmf.y_sealevel
        file["base_level_shift"] = markers_xdmf.base_level_shift

        println(">> RAM before IDs: $(round(Sys.maxrss()/1024/1024, digits=2)) MB")
        let
            group = JLD2.Group(file, jld_dataname_ids)
            group["array"] = collect(Float64, 0:nmarkers-1)
            group["name"] = "marker_ids"
            group["units"] = "None"
        end
        # Reclaim the IDs buffer before the X write burst so the two don't stack.
        GC.gc(false)

        println(">> RAM before X: $(round(Sys.maxrss()/1024/1024, digits=2)) MB")
        let
            marker_x_km = getoutform(markers_xdmf.marker_x_obj)
            marker_x_km ./= 1000.0
            group = JLD2.Group(file, jld_dataname_x)
            group["array"] = marker_x_km
            group["name"] = "marker_x"
            group["units"] = "km"
        end
        # Reclaim the X buffer before the Y write.
        GC.gc(false)

        println(">> RAM before Y: $(round(Sys.maxrss()/1024/1024, digits=2)) MB")
        let
            marker_y_km = getoutform(markers_xdmf.marker_y_obj)
            marker_y_km .*= -1.0 / 1000.0
            group = JLD2.Group(file, jld_dataname_y)
            group["array"] = marker_y_km
            group["name"] = "marker_y"
            group["units"] = "km"
        end
        # Reclaim the Y buffer before the scalar loop.
        GC.gc(false)

        println(">> RAM before scalar loop: $(round(Sys.maxrss()/1024/1024, digits=2)) MB")
        iscalar = 0
        for (obj, meta) in zip(markers_xdmf.enabled_marker_objs, markers_xdmf.scalar_metas)
            group = JLD2.Group(file, meta.jld_dataname)
            group["array"] = getoutform(obj)
            group["name"] = meta.name
            group["units"] = meta.units
            iscalar += 1
            if iscalar % 5 == 0 || iscalar == length(markers_xdmf.scalar_metas)
                println(">> RAM scalar $(iscalar)/$(length(markers_xdmf.scalar_metas)) [$(meta.name)]: $(round(Sys.maxrss()/1024/1024, digits=2)) MB")
            end
        end
        println(">> RAM after scalar loop (pre-close): $(round(Sys.maxrss()/1024/1024, digits=2)) MB")
    end
    println(">> RAM after jldopen close: $(round(Sys.maxrss()/1024/1024, digits=2)) MB")
end

end