module ExtensionInflowAndOutflowAlongSides

import EarthBox.ModelDataContainer: ModelData
import ..Extension: check_required_types_and_domains
import ..Extension: calculate_recycled_air_locations
import ..ResetMarkers: reset_recycled_marker_properties!
import ..ResetMarkers: SubsurfaceProps
import ..ResetMarkers: make_recycle_property_arrays
import ..Coordinates: RecycledCoordinates
import ..Coordinates: MarkerRecycleArrays
import ..Coordinates: InflowOutflowInputs
import ..Coordinates: reset_recycled_marker_coordinates
import ..Coordinates: update_recycled_coordinates_vertical_inflow_outflow_boundary!
import ..InitializeMarkerRecycling: initialize_recycling
import ..CountRecycle: get_left_and_right_recycle_counts
import ...BoundaryVelocity.InflowRamp: get_top_and_bottom_of_left_or_right_ramp
import ...BoundaryVelocity.InflowGeometry: get_starting_inflow_coordinates_left
import ...BoundaryVelocity.InflowGeometry: get_starting_inflow_coordinates_right
import ...BoundaryVelocity.InflowGeometry: get_y_height_of_inflow

function recycle_markers!(model::ModelData, inside_flags::Vector{Int8})::Nothing
    check_required_types_and_domains(model)
    execute_recycling_steps(model, inside_flags)
    return nothing
end

function execute_recycling_steps(model::ModelData, inside_flags::Vector{Int8})::Nothing
    (
        outside_indices, nrecycle_air, nrecycle_subsurface
    ) = initialize_recycling(model, inside_flags)
    recycled_coordinates = calculate_recycled_locations(
        model, nrecycle_air, nrecycle_subsurface)
    marker_recycling = reset_recycled_marker_coordinates(
        model, recycled_coordinates, outside_indices)
    subsurface_props = get_subsurface_reset_properties(model, marker_recycling)
    reset_recycled_marker_properties!(model, marker_recycling, subsurface_props)
    return nothing
end

function calculate_recycled_locations(
    model::ModelData,
    nrecycle_air::Int,
    nrecycle_subsurface::Int
)::RecycledCoordinates
    newx_air, newy_air = calculate_recycled_air_locations(nrecycle_air, model)
    (
        newx_subsurface, newy_subsurface
    ) = calculate_recycled_subsurface_locations(nrecycle_subsurface, model)
    return RecycledCoordinates(newx_air, newy_air, newx_subsurface, newy_subsurface)
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
    newx_subsurface = zeros(nrecycle_subsurface)
    newy_subsurface = zeros(nrecycle_subsurface)

    update_recycled_subsurface_locations_left_side!(
        model, nrecycle_subsurface, use_asymmetric_contraction,
        newx_subsurface, newy_subsurface)

    update_recycled_subsurface_locations_right_side!(
        model, nrecycle_subsurface, use_asymmetric_contraction,
        newx_subsurface, newy_subsurface)

    return newx_subsurface, newy_subsurface
end

function update_recycled_subsurface_locations_left_side!(
    model::ModelData,
    nrecycle_subsurface::Int,
    use_asymmetric_contraction::Bool,
    newx_subsurface::Vector{Float64},
    newy_subsurface::Vector{Float64}
)::Nothing
    nrecycle_left, _ = get_left_and_right_recycle_counts(
        nrecycle_subsurface, use_asymmetric_extension=use_asymmetric_contraction)
    x_initial, y_initial = get_starting_inflow_coordinates_left(model)
    y_smooth_top, y_smooth_bottom = get_top_and_bottom_of_left_or_right_ramp(model, left_side=true)
    inflow_outflow_inputs = InflowOutflowInputs(
        nrecycle        = nrecycle_left, 
        nrecycle_start  = 0, 
        y_height        = get_y_height_of_inflow(model, y_initial),
        y_initial       = y_initial, 
        x_initial       = x_initial, 
        y_smooth_top    = y_smooth_top, 
        y_smooth_bottom = y_smooth_bottom,
        displacement_x  = get_inflow_displacement_left_x(model)
    )
    update_recycled_coordinates_vertical_inflow_outflow_boundary!(
        inflow_outflow_inputs, newy_subsurface, newx_subsurface)
    return nothing
end

function update_recycled_subsurface_locations_right_side!(
    model::ModelData,
    nrecycle_subsurface::Int,
    use_asymmetric_contraction::Bool,
    newx_subsurface::Vector{Float64},
    newy_subsurface::Vector{Float64}
)::Nothing
    nrecycle_left, nrecycle_right = get_left_and_right_recycle_counts(
        nrecycle_subsurface, use_asymmetric_extension=use_asymmetric_contraction)
    x_initial, y_initial = get_starting_inflow_coordinates_right(model)
    y_smooth_top, y_smooth_bottom = get_top_and_bottom_of_left_or_right_ramp(model, right_side=true)
    inflow_outflow_inputs = InflowOutflowInputs(
        nrecycle = nrecycle_right,
        nrecycle_start = nrecycle_left,
        y_height = get_y_height_of_inflow(model, y_initial),
        y_initial = y_initial,
        x_initial = x_initial,
        y_smooth_top = y_smooth_top,
        y_smooth_bottom = y_smooth_bottom,
        displacement_x = get_inflow_displacement_right_x(model)
    )
    update_recycled_coordinates_vertical_inflow_outflow_boundary!(
        inflow_outflow_inputs, newy_subsurface, newx_subsurface)
    return nothing
end

"""
    get_inflow_displacement_left_x(model::ModelData)::Float64

Get left inflow displacements for extension models.

# Returns
- `displacement_left_x`: Displacement along the left side of the model
"""
function get_inflow_displacement_left_x(model::ModelData)::Float64
    velocity = model.bcs.parameters.velocity
    velocity_inflow_left = velocity.velocity_inflow_left.value
    timestep = model.timestep.parameters.main_time_loop.timestep.value
    return abs(velocity_inflow_left) * timestep
end

"""
    get_inflow_displacement_right_x(model::ModelData)::Float64

Get right inflow displacements for extension models.

# Returns
- `displacement_right_x`: Displacement along the right side of the model
"""
function get_inflow_displacement_right_x(model::ModelData)::Float64
    velocity = model.bcs.parameters.velocity
    velocity_inflow_right = velocity.velocity_inflow_right.value
    timestep = model.timestep.parameters.main_time_loop.timestep.value
    return -abs(velocity_inflow_right) * timestep
end

"""
    get_subsurface_reset_properties(
        model::ModelData,
        marker_recycling::MarkerRecycleArrays
    )::SubsurfaceProps

Get subsurface properties for resetting recycled subsurface markers.
"""
function get_subsurface_reset_properties(
    model::ModelData,
    marker_recycling::MarkerRecycleArrays
)::SubsurfaceProps
    recycle_indices = marker_recycling.recycle_indices
    recycle_flags = marker_recycling.recycle_flags
    nrecycle_total = length(recycle_indices)

    (
        matids_recycle, pressure_recycle, temperature_recycle
    ) = make_recycle_property_arrays(nrecycle_total)
    matid_asthenosphere = model.materials.dicts.matid_domains["Asthenosphere"]

    Threads.@threads for irecycle in 1:nrecycle_total
        recycle_flag = recycle_flags[irecycle]
        if recycle_flag == 1
            imarker = recycle_indices[irecycle]
            y_marker = model.markers.arrays.location.marker_y.array[imarker]
            ysize = model.grids.parameters.geometry.ysize.value
            adiabatic_gradient = model.heat_equation.parameters.initial_condition.adiabatic_gradient.value
            temperature_bottom = model.bcs.parameters.temperature.temperature_bottom.value

            pressure_reset = calc_pressure(y_marker)
            temperature_reset = calc_temperature(
                y_marker, ysize, temperature_bottom, adiabatic_gradient)

            matids_recycle[irecycle] = matid_asthenosphere
            pressure_recycle[irecycle] = pressure_reset
            temperature_recycle[irecycle] = temperature_reset
        end
    end

    return SubsurfaceProps(matids_recycle, pressure_recycle, temperature_recycle)
end

"""
    calc_pressure(y_marker::Float64)::Float64

Calculate pressure at the base of the model.
"""
function calc_pressure(y_marker::Float64)::Float64
    return 3300.0 * 9.81 * y_marker
end

"""
    calc_temperature(
        y_marker::Float64,
        ysize::Float64,
        temperature_bottom::Float64,
        adiabatic_gradient::Float64
    )::Float64

Calculate temperature at the base of the model.
"""
function calc_temperature(
    y_marker::Float64,
    ysize::Float64,
    temperature_bottom::Float64,
    adiabatic_gradient::Float64
)::Float64
    dy_km = (ysize - y_marker) / 1000.0
    return temperature_bottom - adiabatic_gradient * dy_km
end

end # module 