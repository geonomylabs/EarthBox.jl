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
    # If true, print the RAM usage before and after each step. Note this adds
    # a lot of output to benchmark output.
    print_ram = false

    jld_marker_filename = markers_xdmf.markers2djld.jld_markerfile
    jld_dataname_x = markers_xdmf.markers2djld.jld_dataname_x
    jld_dataname_y = markers_xdmf.markers2djld.jld_dataname_y
    jld_dataname_ids = markers_xdmf.markers2djld.jld_dataname_ids
    nmarkers = markers_xdmf.markers2djld.nmarkers

    jld_marker_file_path = joinpath(output_dir, jld_marker_filename)
    #println(">> RAM make_jld2_file start: $(round(Sys.maxrss()/1024/1024, digits=2)) MB")
    # Python earthbox markers.py: h5py gzip compression_opts=1 — zlib deflate level 1.
    # Import ZlibCompressor here so CodecZlib is loaded with EarthBox (JLD2's dynamic import
    # can fail when Julia is not started with a project that installs CodecZlib).
    JLD2.jldopen(jld_marker_file_path, "w"; compress=ZlibCompressor(level=1)) do file
        file["time"] = markers_xdmf.model_time
        file["time_units"] = markers_xdmf.time_units
        file["y_sealevel"] = markers_xdmf.y_sealevel
        file["base_level_shift"] = markers_xdmf.base_level_shift

        if print_ram
            println(">> RAM before IDs: $(round(Sys.maxrss()/1024/1024, digits=2)) MB")
        end
        let
            group = JLD2.Group(file, jld_dataname_ids)
            group["array"] = collect(Float64, 0:nmarkers-1)
            group["name"] = "marker_ids"
            group["units"] = "None"
        end
        # Reclaim the IDs buffer before the X write burst so the two don't stack.
        GC.gc(false)

        if print_ram
            println(">> RAM before X: $(round(Sys.maxrss()/1024/1024, digits=2)) MB")
        end
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

        if print_ram
            println(">> RAM before Y: $(round(Sys.maxrss()/1024/1024, digits=2)) MB")
        end
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

        # TODO: This loop is where there is a big spike in RAM usage. Here are some optimization 
        # options from Claude Code that could be explored:
        #
        # Optimization option 1 (this option is currently implemented below but untested)
        #
        # 1. Add GC.gc(false) inside the scalar loop, either every iteration or every 2–3. The IDs/X/Y 
        # blocks prove this pattern keeps the watermark flat; the scalar loop is the only place in this 
        # function missing it. Lowest risk, smallest change.
        #
        # 2. Wrap the loop body in a let ... end, matching the IDs/X/Y blocks. A let-scoped group and 
        # getoutform result become unreachable at the end of the block, giving the GC a clearer signal 
        # than a bare for body where bindings are reused across iterations.
        #
        # Optimization option 2
        # 
        # 1. Extract the per-scalar write into a small function — write_scalar!(file, obj, meta). When 
        # the function returns, every local (the getoutform array, the JLD2.Group, any compression 
        # temporaries held in caller frames) is unreachable by construction. Pair with a GC.gc(false) 
        # after the call, or even just rely on Julia's GC having a cleaner view of lifetime.
        #
        # 2. Reusable scratch buffer. If getoutform(obj) always returns the same size/dtype (one Float64 
        # per marker, I'd guess), adding an in-place variant getoutform!(buf, obj) with a single 
        # pre-allocated buf reused across all 11 scalars eliminates 10 of the 11 allocations outright. 
        # This is the cleanest fix — no GC reliance at all — but requires touching getoutform's API. 
        # Note JLD2 will still allocate its own zlib working buffers; the scratch buffer removes the
        # per-iteration input allocation, not the compressor's state.
        #
        # Optimization option 3
        #
        # 1. Reduce the JLD2 chunk size for scalar datasets. Zlib's working set scales with chunk size; 
        # smaller chunks = smaller compressor footprint per write at the cost of a slightly larger 
        # file and more per-chunk overhead. Worth trying only if 1–4 aren't enough.

        if print_ram
            println(">> RAM before scalar loop: $(round(Sys.maxrss()/1024/1024, digits=2)) MB")
        end
        iscalar = 0
        for (obj, meta) in zip(markers_xdmf.enabled_marker_objs, markers_xdmf.scalar_metas)
            let
                group = JLD2.Group(file, meta.jld_dataname)
                group["array"] = getoutform(obj)
                group["name"] = meta.name
                group["units"] = meta.units
            end
            # Reclaim the scalar buffer + compression temporaries before the next
            # iteration so successive scalar writes do not stack in memory.
            GC.gc(false)
            iscalar += 1
            if print_ram
                if iscalar % 5 == 0 || iscalar == length(markers_xdmf.scalar_metas)
                    println(">> RAM scalar $(iscalar)/$(length(markers_xdmf.scalar_metas)) [$(meta.name)]: $(round(Sys.maxrss()/1024/1024, digits=2)) MB")
                end
            end
        end
        if print_ram
            println(">> RAM after scalar loop (pre-close): $(round(Sys.maxrss()/1024/1024, digits=2)) MB")
        end
    end
    if print_ram
        println(">> RAM after jldopen close: $(round(Sys.maxrss()/1024/1024, digits=2)) MB")
    end
end

end