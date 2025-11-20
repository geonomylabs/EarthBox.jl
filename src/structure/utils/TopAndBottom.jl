module TopAndBottom

import EarthBox.ModelDataContainer: ModelData
import ..SmoothSurface: smooth_surface

function calculate_search_radius(
    mxstep::Float64,
    gridx::Vector{Float64},
    marker_search_factor::Float64=2.0
)::Float64
    delta_x_grid = (gridx[2] - gridx[1])/2.0
    search_radius_minimum = mxstep*marker_search_factor
    search_radius = max(delta_x_grid, search_radius_minimum)
    return search_radius
end

function calculate_layer_thickness(
    model::ModelData,
    material_ids_of_layer::Vector{Int16},
    gridx::Vector{Float64};
    use_smoothing::Bool=true
)::Vector{Float64}
    marker_arrays = model.markers.arrays
    marker_matids = marker_arrays.material.marker_matid.array
    location = marker_arrays.location
    marker_x = location.marker_x.array
    marker_y = location.marker_y.array

    marker_search_factor = model.topography.parameters.topo_grid.marker_search_factor.value
    nsmooth_top_bottom = model.topography.parameters.topo_grid.nsmooth_top_bottom.value
    mxstep = model.markers.parameters.distribution.mxstep.value
    search_radius = calculate_search_radius(mxstep, gridx, marker_search_factor)
    top, bottom = calculate_top_and_bottom_of_layer_opt(
        material_ids_of_layer, marker_matids,
        marker_x, marker_y, gridx, search_radius, use_smoothing=false
    )
    if use_smoothing
        # Smooth the sediment top and bottom surfaces to avoid high frequency
        # oscillations in the sediment thickness
        top = smooth_surface(top, nsmooth=nsmooth_top_bottom)
        bottom = smooth_surface(bottom, nsmooth=nsmooth_top_bottom)
    end

    thickness = bottom .- top
    return thickness
end

""" Find the shallowest and deepest y-coordinates of a material layer.

The material layer is a collection of material ids (material_ids).

Note that this function does not take into account vertically discontinuous layers.
"""
function calculate_top_and_bottom_of_layer(
    material_ids_of_layer::Vector{Int16},
    marker_matid::Vector{Int16},
    marker_x::Vector{Float64},
    marker_y::Vector{Float64},
    gridx::Vector{Float64},
    search_radius::Float64
)::Tuple{Vector{Float64}, Vector{Float64}}
    marknum = length(marker_x)
    xnum = length(gridx)
    tops = zeros(Float64, xnum)
    bottoms = zeros(Float64, xnum)
    
    Threads.@threads for j in 1:xnum
        ymin = 1e32
        ymax = -1e32
        xgrid = gridx[j]
        xmin_search = xgrid - search_radius
        xmax_search = xgrid + search_radius
        ifind_top = 0
        ifind_bottom = 0
        
        for imarker in 1:marknum
            xmarker = marker_x[imarker]
            ymarker = marker_y[imarker]
            mid = marker_matid[imarker]
            if xmin_search <= xmarker <= xmax_search
                if mid in material_ids_of_layer && ymarker > ymax
                    ymax = ymarker
                    ifind_bottom = 1
                end
                if mid in material_ids_of_layer && ymarker < ymin
                    ymin = ymarker
                    ifind_top = 1
                end
            end
        end
        
        if ifind_top == 1
            tops[j] = ymin
        end
        if ifind_bottom == 1
            bottoms[j] = ymax
        end
    end
    
    return tops, bottoms
end

""" Find the shallowest and deepest y-coordinates of a material layer.

This is the optimized version of this function.

The material layer is a collection of material ids (material_ids).

Note that this function does not take into account vertically discontinuous layers.
"""
function calculate_top_and_bottom_of_layer_opt(
    material_ids_of_layer::Vector{Int16},
    marker_matid::Vector{Int16},
    marker_x::Vector{Float64},
    marker_y::Vector{Float64},
    gridx::Vector{Float64},
    search_radius::Float64;
    use_smoothing::Bool=true,
    nsmooth::Int=2
)::Tuple{Vector{Float64}, Vector{Float64}}
    marker_indices_layer = filter_markers_outside_of_layer(marker_matid, material_ids_of_layer)
    marknum_layer = length(marker_indices_layer)

    xnum = length(gridx)
    tops = zeros(Float64, xnum)
    bottoms = zeros(Float64, xnum)
    
    Threads.@threads for j in 1:xnum
        ymin = 1e32
        ymax = -1e32
        xgrid = gridx[j]
        xmin_search = xgrid - search_radius
        xmax_search = xgrid + search_radius
        ifind_top = 0
        ifind_bottom = 0
        
        for imarker_layer in 1:marknum_layer
            imarker = marker_indices_layer[imarker_layer]
            xmarker = marker_x[imarker]
            ymarker = marker_y[imarker]
            mid = marker_matid[imarker]
            if xmin_search <= xmarker <= xmax_search
                if mid in material_ids_of_layer && ymarker > ymax
                    ymax = ymarker
                    ifind_bottom = 1
                end
                if mid in material_ids_of_layer && ymarker < ymin
                    ymin = ymarker
                    ifind_top = 1
                end
            end
        end
        
        if ifind_top == 1
            tops[j] = ymin
        end
        if ifind_bottom == 1
            bottoms[j] = ymax
        end
    end

    if use_smoothing
        tops = smooth_surface(tops, nsmooth=nsmooth)
        bottoms = smooth_surface(bottoms, nsmooth=nsmooth)
    end

    return tops, bottoms
end

function filter_markers_outside_of_layer(
    marker_matid::Vector{Int16},
    material_ids_of_layer::Vector{Int16}
)::Vector{Int64}
    marknum = length(marker_matid)
    marker_indices_layer_tmp = Vector{Int64}(undef, marknum)

    icount = 0
    for imarker in 1:marknum
        mid = marker_matid[imarker]
        if mid in material_ids_of_layer
            marker_indices_layer_tmp[icount+1] = imarker
            icount += 1
        else
            marker_indices_layer_tmp[imarker] = -1
        end
    end

    marker_indices_layer = Vector{Int64}(undef, icount)
    
    Threads.@threads for i in 1:icount
        marker_indices_layer[i] = marker_indices_layer_tmp[i]
    end

    return marker_indices_layer
end

""" Find the shallowest and deepest y-coordinates of a material layer.

The material layer is a collection of material ids (material_ids).

Note that this function does not take into account vertically discontinuous layers.

Optimization version of the function below.
"""
function calculate_top_and_bottom_of_swarm_opt(
    x_sorted_marker_indices_swarm::Vector{Int64},
    marker_x::Vector{Float64},
    marker_y::Vector{Float64},
    gridx::Vector{Float64},
    search_radius::Float64;
    use_smoothing::Bool=true,
    nsmooth::Int=2
)::Tuple{Vector{Float64}, Vector{Float64}}
    nswarm = length(x_sorted_marker_indices_swarm)
    xnum = length(gridx)
    tops = fill(1e32, xnum)
    bottoms = fill(-1e32, xnum)
    
    Threads.@threads for m in 1:nswarm
        imarker = x_sorted_marker_indices_swarm[m]
        xmarker = marker_x[imarker]
        ymarker = marker_y[imarker]
        # Binary search for the range of grid points affected by this marker
        left_gridx_index = searchsortedfirst(gridx, xmarker - search_radius)
        right_gridx_index = searchsortedlast(gridx, xmarker + search_radius)
        for j in left_gridx_index:right_gridx_index
            if abs(xmarker - gridx[j]) <= search_radius
                if ymarker > bottoms[j]
                    bottoms[j] = ymarker
                end
                if ymarker < tops[j]
                    tops[j] = ymarker
                end
            end
        end
    end

    Threads.@threads for j in 1:xnum
        if tops[j] == 1e32
            tops[j] = 0.0
        end
        if bottoms[j] == -1e32
            bottoms[j] = 0.0
        end
    end

    if use_smoothing
        tops = smooth_surface(tops, nsmooth=nsmooth)
        bottoms = smooth_surface(bottoms, nsmooth=nsmooth)
    end

    return tops, bottoms
end

""" Find the shallowest and deepest y-coordinates of a material layer.

The material layer is a collection of material ids (material_ids).

Note that this function does not take into account vertically discontinuous layers.
"""
function calculate_top_and_bottom_of_swarm(
    marker_indices_swarm::Vector{Int64},
    marker_x::Vector{Float64},
    marker_y::Vector{Float64},
    gridx::Vector{Float64},
    search_radius::Float64;
    use_smoothing::Bool=true,
    nsmooth::Int=2
)::Tuple{Vector{Float64}, Vector{Float64}}
    nswarm = length(marker_indices_swarm)
    xnum = length(gridx)
    tops = zeros(Float64, xnum)
    bottoms = zeros(Float64, xnum)
    
    Threads.@threads for j in 1:xnum
        ymin = 1e32
        ymax = -1e32
        xgrid = gridx[j]
        xmin_search = xgrid - search_radius
        xmax_search = xgrid + search_radius
        ifind_top = 0
        ifind_bottom = 0
        
        for m in 1:nswarm
            imarker = marker_indices_swarm[m]
            xmarker = marker_x[imarker]
            ymarker = marker_y[imarker]
            if xmin_search <= xmarker <= xmax_search
                if ymarker > ymax
                    ymax = ymarker
                    ifind_bottom = 1
                end
                if ymarker < ymin
                    ymin = ymarker
                    ifind_top = 1
                end
            end
        end
        
        if ifind_top == 1
            tops[j] = ymin
        end
        if ifind_bottom == 1
            bottoms[j] = ymax
        end
    end

    if use_smoothing
        tops = smooth_surface(tops, nsmooth=nsmooth)
        bottoms = smooth_surface(bottoms, nsmooth=nsmooth)
    end

    return tops, bottoms
end

end