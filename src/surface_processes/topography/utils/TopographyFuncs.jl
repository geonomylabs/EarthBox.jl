module TopographyFuncs

import EarthBox.ModelDataContainer: ModelData
import EarthBox.MathTools: linear_interp_at_x_location, linear_interp_bisection
import EarthBox.Markers.MarkerMaterials.InitManager: FlexureTriangularHole

function get_topo_coordinates(
    gridt::Matrix{Float64}
)::Tuple{Vector{Float64},Vector{Float64}}
    toponum = size(gridt, 2)
    topo_gridx = zeros(toponum)
    topo_gridy = zeros(toponum)
    for j in 1:toponum
        topo_gridx[j] = gridt[1, j]
        topo_gridy[j] = gridt[2, j]
    end
    return topo_gridx, topo_gridy
end

function topography_and_carb_grid_initialize!(model::ModelData)::Nothing
    initialize_topo_grid!(model, model.topography.arrays.gridt.array)
    initialize_topo_grid!(model, model.carbonate.arrays.grid_carb.array)
    return nothing
end

function topography_and_carb_grid_initialize_triangular_hole!(
    model::ModelData,
    use_crustal_hole::Bool
)::Nothing
    if use_crustal_hole
        initialize_topo_array_for_crustal_hole(
            model, model.topography.arrays.gridt.array)
        initialize_topo_array_for_crustal_hole(
            model, model.carbonate.arrays.grid_carb.array)
    end
    return nothing
end

function initialize_topo_grid!(model::ModelData, gridt::Matrix{Float64})::Nothing
    topo_grid = model.topography.parameters.topo_grid
    toponum = topo_grid.toponum.value
    dx_topo = topo_grid.dx_topo.value
    depth_sticky = model.geometry.parameters.sticky_air_geometry.thick_air.value

    # Horizontal grid for topography
    gridt[1, 1] = 0.0  # beginning of topography profile
    for i in 2:toponum
        gridt[1, i] = gridt[1, i-1] + dx_topo
    end

    # Initial elevation for topography profile
    for i in 1:toponum
        gridt[2, i] = depth_sticky
    end
    return nothing
end

function initialize_topo_array_for_crustal_hole(
    model::ModelData,
    gridt::Matrix{Float64}
)::Nothing
    toponum = model.topography.parameters.topo_grid.toponum.value
    depth_sticky = model.geometry.parameters.sticky_air_geometry.thick_air.value
    xhole_start = model.geometry.parameters.crustal_hole.xhole_start.value
    xhole_middle = model.geometry.parameters.crustal_hole.xhole_middle.value
    xhole_end = model.geometry.parameters.crustal_hole.xhole_end.value
    xhole_depth = model.geometry.parameters.crustal_hole.xhole_depth.value
    # Initial elevation for topography profile
    for i in 1:toponum
        adjust_topography_for_hole(
            i, gridt, depth_sticky,
            xhole_start, xhole_middle, xhole_end, xhole_depth
        )
    end
    return nothing
end

function adjust_topography_for_hole(
    i::Int,
    gridt::Matrix{Float64},
    thick_water::Float64,
    xhole_start::Float64,
    xhole_middle::Float64,
    xhole_end::Float64,
    xhole_depth::Float64
)::Nothing
    marker_x = gridt[1, i]
    marker_y = thick_water

    on_left_side = FlexureTriangularHole.on_left_side_of_hole(
        marker_x, marker_y, xhole_start, xhole_middle, thick_water)
    on_right_side = FlexureTriangularHole.on_right_side_of_hole(
        marker_x, marker_y, xhole_middle, xhole_end, thick_water)

    if on_left_side
        gridt[2, i] = FlexureTriangularHole.calc_depth_hole_left(
            marker_x, thick_water, xhole_depth, xhole_start, xhole_middle)
    end
    if on_right_side
        gridt[2, i] = FlexureTriangularHole.calc_depth_hole_right(
            marker_x, thick_water, xhole_depth, xhole_middle, xhole_end)
    end
    return nothing
end

function get_carbonate_growth_depth(
    y_depths::Vector{Float64},
    depth_sticky::Float64,
    z_photic_m::Float64
)::Tuple{Float64, Float64, Float64}
    # Find minima and maxima indices
    diff_sign = diff(sign(diff(y_depths)))
    imin_all = findall(diff_sign .> 0) .+ 1
    imax_all = findall(diff_sign .< 0) .+ 1
    
    # Find the maximum depth of the maxima within photic zone
    ymax = -1e32
    for i in eachindex(imax_all)
        ycheck = y_depths[imax_all[i]]
        if z_photic_m >= ycheck >= depth_sticky && ycheck > ymax
            ymax = ycheck
        end
    end
    
    # Find minimum depth of the minima within photic zone
    ymin = depth_sticky
    for i in eachindex(imin_all)
        ycheck = y_depths[imin_all[i]]
        if z_photic_m >= ycheck >= depth_sticky && ycheck < ymin
            ymin = ycheck
        end
    end
    
    # Calculate carbonate growth depth
    zmax_carb = depth_sticky
    if ymax < depth_sticky && ymin > z_photic_m
        zmax_carb = depth_sticky
    elseif ymax < depth_sticky && ymin < z_photic_m
        ymax = z_photic_m
        a_factor = (ymax - ymin) / 2
        zmax_carb = ymin + a_factor
    elseif ymax > depth_sticky && ymin < z_photic_m
        a_factor = (ymax - ymin) / 2
        zmax_carb = ymin + a_factor
    end
    
    return ymax, ymin, zmax_carb
end

end # module TopographyFuncs 