module ContractionRecycleLocations

import EarthBox.ModelDataContainer: ModelData
import ...BoundaryVelocity.StickyThickness: get_sticky_thickness
import ..Coordinates: RecycledCoordinates
import ..Coordinates: InflowOutflowInputs
import ..Coordinates: update_recycled_coordinates_vertical_boundary!
import ..CountRecycle: get_left_and_right_recycle_counts

function calculate_recycled_locations(
    model::ModelData,
    nrecycle_air::Int,
    nrecycle_subsurface::Int,
    use_asymmetric_contraction::Bool = false
)::RecycledCoordinates
    newx_air, newy_air = calculate_recycled_air_locations(
        nrecycle_air, model, use_asymmetric_contraction)
    newx_subsurface, newy_subsurface = calculate_recycled_subsurface_locations(
        nrecycle_subsurface, model, use_asymmetric_contraction)
    return RecycledCoordinates(newx_air, newy_air, newx_subsurface, newy_subsurface)
end

"""
    calculate_recycled_air_locations(
        nrecycle_air::Int,
        model::ModelData,
        use_asymmetric_contraction::Bool = false
    )::Tuple{Vector{Float64}, Vector{Float64}}

Calculate recycled air locations.

# Returns
- `newx_air`: New recycled x-locations for sticky air/water markers
- `newy_air`: New recycled y-locations for sticky air/water markers
"""
function calculate_recycled_air_locations(
    nrecycle_air::Int,
    model::ModelData,
    use_asymmetric_contraction::Bool = false
)::Tuple{Vector{Float64}, Vector{Float64}}
    displacement_x = get_displacement_velocity(model, use_asymmetric_contraction)

    nrecycle_left, nrecycle_right = get_left_and_right_recycle_counts(
        nrecycle_air, use_asymmetric_extension=use_asymmetric_contraction)

    (
        sticky_thickness_left, sticky_thickness_right
    ) = get_sticky_thickness(model)

    newx_air = zeros(nrecycle_air)
    newy_air = zeros(nrecycle_air)

    # Left Boundary
    inflow_outflow_inputs = InflowOutflowInputs(
        nrecycle = nrecycle_left,
        nrecycle_start = 0,
        y_height = sticky_thickness_left,
        y_initial = 0.0,
        x_initial = 0.0,
        y_smooth_top = 0.0,
        y_smooth_bottom = 0.0,
        displacement_x = displacement_x
    )

    update_recycled_coordinates_vertical_boundary!(
        inflow_outflow_inputs, newy_air, newx_air)

    # Right Boundary
    xsize = model.grids.parameters.geometry.xsize.value
    inflow_outflow_inputs = InflowOutflowInputs(
        nrecycle = nrecycle_right,
        nrecycle_start = nrecycle_left,
        y_height = sticky_thickness_right,
        y_initial = 0.0,
        x_initial = xsize,
        y_smooth_top = 0.0,
        y_smooth_bottom = 0.0,
        displacement_x = -displacement_x
    )

    update_recycled_coordinates_vertical_boundary!(
        inflow_outflow_inputs, newy_air, newx_air)

    return newx_air, newy_air
end

"""
    calculate_recycled_subsurface_locations(
        nrecycle_subsurface::Int,
        model::ModelData,
        use_asymmetric_contraction::Bool = false
    )::Tuple{Vector{Float64}, Vector{Float64}}

Calculate recycled lithosphere locations.

# Returns
- `newx_subsurface`: New recycled x-locations for subsurface markers
- `newy_subsurface`: New recycled y-locations for subsurface markers
"""
function calculate_recycled_subsurface_locations(
    nrecycle_subsurface::Int,
    model::ModelData,
    use_asymmetric_contraction::Bool = false
)::Tuple{Vector{Float64}, Vector{Float64}}
    displacement_x = get_displacement_velocity(model, use_asymmetric_contraction)
    nrecycle_left, nrecycle_right = get_left_and_right_recycle_counts(
        nrecycle_subsurface, use_asymmetric_extension=use_asymmetric_contraction)

    (
        sticky_thickness_left, sticky_thickness_right
    ) = get_sticky_thickness(model)

    newx_subsurface = zeros(nrecycle_subsurface)
    newy_subsurface = zeros(nrecycle_subsurface)

    # Left Boundary
    ysize = model.grids.parameters.geometry.ysize.value
    inflow_outflow_inputs = InflowOutflowInputs(
        nrecycle = nrecycle_left,
        nrecycle_start = 0,
        y_height = ysize - sticky_thickness_left,
        y_initial = sticky_thickness_left,
        x_initial = 0.0,
        y_smooth_top = 0.0,
        y_smooth_bottom = 0.0,
        displacement_x = displacement_x
    )

    update_recycled_coordinates_vertical_boundary!(
        inflow_outflow_inputs, newy_subsurface, newx_subsurface)

    # Right Boundary
    inflow_outflow_inputs = InflowOutflowInputs(
        nrecycle = nrecycle_right,
        nrecycle_start = nrecycle_left,
        y_height = ysize - sticky_thickness_right,
        y_initial = sticky_thickness_right,
        x_initial = model.grids.parameters.geometry.xsize.value,
        y_smooth_top = 0.0,
        y_smooth_bottom = 0.0,
        displacement_x = -displacement_x
    )

    update_recycled_coordinates_vertical_boundary!(
        inflow_outflow_inputs, newy_subsurface, newx_subsurface)

    return newx_subsurface, newy_subsurface
end

"""
    get_displacement_velocity(
        model::ModelData,
        use_asymmetric_contraction::Bool = false
    )::Float64

Get displacement velocity for contraction models.

# Returns
- `displacement_side`: Displacement velocity for contraction models
"""
function get_displacement_velocity(
    model::ModelData,
    use_asymmetric_contraction::Bool = false
)::Float64
    timestep = model.timestep.parameters.main_time_loop.timestep.value
    full_velocity_contraction = model.bcs.parameters.velocity.full_velocity_contraction.value

    if use_asymmetric_contraction
        displacement_side = abs(full_velocity_contraction) * timestep
    else
        displacement_side = abs(full_velocity_contraction) * timestep * 0.5
    end
    return displacement_side
end

end # module 