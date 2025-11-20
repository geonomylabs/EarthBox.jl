module Extension

import EarthBox.ModelDataContainer: ModelData
import ..ResetMarkers: reset_recycled_marker_properties!
import ..ResetMarkers: make_recycle_property_arrays
import ..ResetMarkers: SubsurfaceProps
import ..Coordinates: reset_recycled_marker_coordinates
import ..Coordinates: calc_recycled_x_coordinates_horizontal_boundary
import ..Coordinates: calc_recycled_y_coordinates_horizontal_boundary
import ..Coordinates: RecycledCoordinates
import ..Coordinates: MarkerRecycleArrays
import ..CheckDicts: check_domain_dict
import ..PressureReset: calculate_reset_pressure
import ..InitializeMarkerRecycling: initialize_recycling

function recycle_markers!(model::ModelData, inside_flags::Vector{Int8})::Nothing
    check_required_types_and_domains(model)
    execute_recycling_steps(model, inside_flags)
    return nothing
end

function check_required_types_and_domains(model::ModelData)::Nothing
    matid_types = model.materials.dicts.matid_types
    type_dict = Dict(
        "StickyAir" => matid_types["StickyAir"][1],
        "StickyWater" => matid_types["StickyWater"][1]
    )
    check_domain_dict(type_dict, "type", "extensional")

    domains = model.materials.dicts.matid_domains
    domain_dict = Dict(
        "Asthenosphere" => domains["Asthenosphere"]
    )
    check_domain_dict(domain_dict, "domain", "extensional")
    return nothing
end

function execute_recycling_steps(model::ModelData, inside_flags::Vector{Int8})::Nothing
    outside_indices, nrecycle_air, nrecycle_subsurface = initialize_recycling(model, inside_flags)
    recycled_coordinates = calculate_recycled_locations(model, nrecycle_air, nrecycle_subsurface)
    marker_recycling = reset_recycled_marker_coordinates(model, recycled_coordinates, outside_indices)
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
    newx_subsurface, newy_subsurface = calculate_recycled_subsurface_locations(
        nrecycle_subsurface, model)
    return RecycledCoordinates(newx_air, newy_air, newx_subsurface, newy_subsurface)
end

function calculate_recycled_air_locations(
    nrecycle_air::Int,
    model::ModelData
)::Tuple{Vector{Float64}, Vector{Float64}}
    newy_air = calc_recycled_y_coordinates_horizontal_boundary(
        nrecycle_air, 0.0, calc_displacement_y_upper(model))
    newx_air = calc_recycled_x_coordinates_horizontal_boundary(
        nrecycle_air, 0.0, model.grids.parameters.geometry.xsize.value)
    return newx_air, newy_air
end

function calc_displacement_y_upper(model::ModelData)::Float64
    velocity_y_upper = model.bcs.parameters.velocity.vyu.value
    timestep = model.timestep.parameters.main_time_loop.timestep.value
    return abs(velocity_y_upper * timestep)
end

function calculate_recycled_subsurface_locations(
    nrecycle_subsurface::Int,
    model::ModelData
)::Tuple{Vector{Float64}, Vector{Float64}}
    xsize = model.grids.parameters.geometry.xsize.value
    ysize = model.grids.parameters.geometry.ysize.value
    newx_subsurface = calc_recycled_x_coordinates_horizontal_boundary(
        nrecycle_subsurface, 0.0, xsize)
    displacement_lower = calc_displacement_y_lower(model)
    newy_subsurface = calc_recycled_y_coordinates_horizontal_boundary(
        nrecycle_subsurface, ysize, displacement_lower)
    return newx_subsurface, newy_subsurface
end

function calc_displacement_y_lower(model::ModelData)::Float64
    velocity_y_lower = model.bcs.parameters.velocity.vyl.value
    timestep = model.timestep.parameters.main_time_loop.timestep.value
    return -abs(velocity_y_lower * timestep)
end

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
            pressure_reset = calc_pressure(model)
            temperature_reset = calc_temperature(model)
            matids_recycle[irecycle] = matid_asthenosphere
            pressure_recycle[irecycle] = pressure_reset
            temperature_recycle[irecycle] = temperature_reset
        end
    end

    return SubsurfaceProps(matids_recycle, pressure_recycle, temperature_recycle)
end

function calc_pressure(model::ModelData)::Float64
    ynum = model.grids.parameters.geometry.ynum.value
    gridy_max = model.grids.arrays.basic.gridy_b.array[ynum]
    return calculate_reset_pressure(gridy_max)
end

function calc_temperature(model::ModelData)::Float64
    return model.bcs.parameters.temperature.temperature_bottom.value
end

end # module Extension 