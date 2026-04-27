module GridFuncs

import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer.Grids2dContainer.ParameterCollection.GridGeometryGroup: GridGeometry

function get_grid_info(
    model::ModelData,
    grid_type::String
)::Tuple{Vector{Float64}, Vector{Float64}, String}
    data = model.grids.arrays
    grid_type = model.grids.parameters.geometry.grid_type
    if grid_type == "basic"
        gridx = data.basic.gridx_b.array
        gridy = data.basic.gridy_b.array
        length_units = data.basic.gridx_b.units
    elseif grid_type == "pressure"
        gridx = data.pressure.gridx_pr.array
        gridy = data.pressure.gridy_pr.array
        length_units = data.pressure.gridx_pr.units
    else
        throw(ArgumentError("grid_type $grid_type is not valid."))
    end
    return (gridx, gridy, length_units)
end

function get_indices_of_markers_outside_domain(
    model::ModelData,
    inside_flags::Vector{Int8}
)::Vector{Int64}
    marknum = model.markers.parameters.distribution.marknum.value
    # Use persistent marknum-sized scratch instead of allocating a fresh
    # Vector{Int64}(marknum) per call. The threaded fill marks each slot
    # with `i` (outside) or `-1` (inside) in marker-index order. The
    # subsequent single-pass packed-prefix loop compacts the outside
    # entries to the front of the same buffer; since `nrecycle <= i`
    # throughout, in-place writes never clobber unread positions. The
    # returned tight Vector preserves the original marker-index order
    # (matching the previous `findall + index` semantics).
    scratch = model.markers.arrays.recycle.marker_outside_indices_scratch.array
    Threads.@threads for i in 1:marknum
        scratch[i] = inside_flags[i] != 1 ? i : -1
    end
    nrecycle = 0
    @inbounds for i in 1:marknum
        v = scratch[i]
        if v >= 0
            nrecycle += 1
            scratch[nrecycle] = v
        end
    end
    return scratch[1:nrecycle]
end

function get_marker_inside_flags(model::ModelData)::Vector{Int8}
    marknum = model.markers.parameters.distribution.marknum.value
    location = model.markers.arrays.location
    marker_x = location.marker_x.array
    marker_y = location.marker_y.array
    # Use the persistent marknum-sized scratch buffer instead of allocating
    # a fresh Vector{Int8}(marknum) per call. The buffer is always fully
    # overwritten below (every index 1..marknum gets either 1 or -1), so
    # stale values from a prior call do not leak through. All downstream
    # consumers in the time loop treat the result as read-only; the
    # threaded fill below is the only writer. The two call sites in
    # TimeLoop.run_loop! / execute_post_solver_steps! are sequential, so
    # there is no race on the shared buffer.
    inside_flags = model.markers.arrays.recycle.marker_inside_flags_buffer.array
    xmin, xmax, ymin, ymax = get_domain_limits(model.grids.parameters.geometry)
    Threads.@threads for i in 1:marknum
        x_marker = marker_x[i]
        y_marker = marker_y[i]
        in_domain = check_in_domain_optimized(xmin, xmax, ymin, ymax, x_marker, y_marker)
        if in_domain
            inside_flags[i] = 1
        else
            inside_flags[i] = -1
        end
    end
    return inside_flags
end

function check_in_domain(
    geometry::GridGeometry,
    x_marker::Float64,
    y_marker::Float64
)::Bool
    gx_min = geometry.xmin.value
    gx_max = geometry.xmax.value
    gy_min = geometry.ymin.value
    gy_max = geometry.ymax.value
    return gx_min <= x_marker <= gx_max && gy_min <= y_marker <= gy_max
end

function get_domain_limits(geometry::GridGeometry)
    gx_min = geometry.xmin.value
    gx_max = geometry.xmax.value
    gy_min = geometry.ymin.value
    gy_max = geometry.ymax.value
    return (gx_min, gx_max, gy_min, gy_max)
end

function check_in_domain_optimized(
    gx_min::Float64,
    gx_max::Float64,
    gy_min::Float64,
    gy_max::Float64,
    x_marker::Float64,
    y_marker::Float64
)
    return gx_min <= x_marker <= gx_max && gy_min <= y_marker <= gy_max
end

end # module 