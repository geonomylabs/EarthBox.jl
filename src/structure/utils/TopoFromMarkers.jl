module TopoFromMarkers

import EarthBox.ModelDataContainer: ModelData
import ..SmoothSurface: smooth_surface
import ..TopAndBottom

""" Calculate topography from markers.

# Updated Arrays
- `model.topography.arrays.gridt`: Array{Float64,2}(7, toponum)
    The sub-array gridt[2,toponum], which stores the y-coordinate of topography is updated.
"""
function calc_topo_from_markers!(model::ModelData)
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marknum = model.markers.parameters.distribution.marknum.value
    marker_matid = model.markers.arrays.material.marker_matid.array

    gridt = model.topography.arrays.gridt.array
    toponum = model.topography.parameters.topo_grid.toponum.value
    topo_gridx = copy(gridt[1,:])

    marker_search_factor = model.topography.parameters.topo_grid.marker_search_factor.value
    mxstep = model.markers.parameters.distribution.mxstep.value
    search_radius = TopAndBottom.calculate_search_radius(
        mxstep, topo_gridx, marker_search_factor)

    matids_sticky_air = model.materials.dicts.matid_types["StickyAir"]
    matids_sticky_water = model.materials.dicts.matid_types["StickyWater"]

    matid_air = matids_sticky_air[1]
    matid_water = matids_sticky_water[1]

    matids_sticky = zeros(Int16, 2)
    matids_sticky[1] = matid_air
    matids_sticky[2] = matid_water

    # Define a dummy array that does not include sticky ids
    matids_rocks = fill(-1, 200)
    for i in 1:100
        if !(i in matids_sticky)
            matids_rocks[i] = i
        end
    end

    thick_air = model.geometry.parameters.sticky_air_geometry.thick_air.value
    max_bathymetry = 25000.0

    depth_max = thick_air + max_bathymetry
    marknum_new, marker_x_new, marker_y_new, marker_matid_new = 
        filter_markers_based_on_depth(
            depth_max, marknum, marker_x, marker_y, marker_matid)

    nsmooth_top_bottom = model.topography.parameters.topo_grid.nsmooth_top_bottom.value

    _, ybottoms_sticky_water = TopAndBottom.calculate_top_and_bottom_of_layer_opt(
        matids_sticky, marker_matid_new,
        marker_x_new, marker_y_new,
        topo_gridx, search_radius,
        use_smoothing=true, nsmooth=nsmooth_top_bottom
    )

    ytops_rock, _ = TopAndBottom.calculate_top_and_bottom_of_layer_opt(
        matids_rocks, marker_matid_new,
        marker_x_new, marker_y_new,
        topo_gridx, search_radius,
        use_smoothing=true, nsmooth=nsmooth_top_bottom
    )

    ytopo = calc_topography(toponum, ytops_rock, ybottoms_sticky_water)

    topo_gridy = smooth_surface(ytopo, nsmooth=nsmooth_top_bottom)
    for i in 1:toponum
        gridt[2,i] = topo_gridy[i]
    end
end

""" Filter markers based on depth.

# Returns
- `marknum_new`: Number of markers after filtering
- `marker_x_new`: Filtered x-coordinates of markers
- `marker_y_new`: Filtered y-coordinates of markers
- `marker_matid_new`: Filtered material IDs of markers
"""
function filter_markers_based_on_depth(
    depth_max::Float64,
    marknum::Int64,
    marker_x::Vector{Float64},
    marker_y::Vector{Float64},
    marker_matid::Vector{Int16}
)
    marker_x_tmp = ones(Float64, marknum)
    marker_y_tmp = ones(Float64, marknum)
    marker_matid_tmp = ones(Int16, marknum)
    marknum_new = 0
    
    for imarker in 1:marknum
        y_marker = marker_y[imarker]
        if y_marker <= depth_max
            x_marker = marker_x[imarker]
            matid = marker_matid[imarker]
            marknum_new += 1
            marker_x_tmp[marknum_new] = x_marker
            marker_y_tmp[marknum_new] = y_marker
            marker_matid_tmp[marknum_new] = matid
        end
    end
    
    marker_x_new = ones(Float64, marknum_new)
    marker_y_new = ones(Float64, marknum_new)
    marker_matid_new = ones(Int16, marknum_new)
    
    for imarker in 1:marknum_new
        marker_x_new[imarker] = marker_x_tmp[imarker]
        marker_y_new[imarker] = marker_y_tmp[imarker]
        marker_matid_new[imarker] = marker_matid_tmp[imarker]
    end

    return marknum_new, marker_x_new, marker_y_new, marker_matid_new
end

""" Calculate topography.

# Arguments
- `xnum`: Number of x-coordinates in topography grid
- `ytops_rock`: Y-coordinates of top of rock layer
- `ybottoms_sticky_water`: Y-coordinates of bottom of sticky water layer

# Returns
- `ytopo`: Y-coordinates of topography (average of top of rock and bottom of sticky water)
"""
function calc_topography(
    xnum::Int64,
    ytops_rock::Vector{Float64},
    ybottoms_sticky_water::Vector{Float64}
)
    ytopo = zeros(Float64, xnum)
    for i in 1:xnum
        ytopo[i] = (ytops_rock[i] + ybottoms_sticky_water[i]) / 2.0
    end
    return ytopo
end

end # module 