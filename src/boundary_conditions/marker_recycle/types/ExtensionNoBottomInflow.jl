module ExtensionNoBottomInflow

import EarthBox.ModelDataContainer: ModelData
import EarthBox: GridFuncs
import ..CheckDicts: check_domain_dict
import ..ResetMarkers: reset_marker!

function recycle_markers!(model::ModelData, inside_flags::Vector{Int8})::Nothing
    check_required_types_and_domains(model)
    execute_recycling_steps(model, inside_flags)
    return nothing
end

function check_required_types_and_domains(model::ModelData)::Nothing
    matid_types = model.materials.dicts.matid_types
    type_dict = Dict(
        "StickyAir" => matid_types["StickyAir"][1]
    )
    check_domain_dict(type_dict, "type", "extensional_no_bottom_flow")
    return nothing
end

function execute_recycling_steps(model::ModelData, inside_flags::Vector{Int8})::Nothing
    outside_indices = GridFuncs.get_indices_of_markers_outside_domain(model, inside_flags)
    reset_recycled_marker_coordinates!(model, outside_indices)
    reset_recycled_marker_properties!(model, outside_indices)
    return nothing
end

"""
    reset_recycled_marker_coordinates!(
        model::ModelData,
        outside_indices::Vector{Int64}
    )::Nothing

Reset recycled marker coordinates.

# Arguments
- `model`: Model data container
- `outside_indices`: Indices of markers outside the model domain

# Updated arrays from group `model.markers.arrays.location`
- `marker_x.array::Vector{Float64}`: X-locations of markers in meters.
- `marker_y.array::Vector{Float64}`: Y-locations of markers in meters.
"""
function reset_recycled_marker_coordinates!(
    model::ModelData,
    outside_indices::Vector{Int64}
)::Nothing
    nrecycle = length(outside_indices)
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array

    newx, newy = calculate_recycled_locations_along_top_boundary(nrecycle, model)

    for ioutside in 1:nrecycle
        imarker = outside_indices[ioutside]
        marker_x[imarker] = newx[ioutside]
        marker_y[imarker] = newy[ioutside]
    end
    return nothing
end

"""
    reset_recycled_marker_properties!(
        model::ModelData,
        outside_indices::Vector{Int64}
    )::Nothing

Reset recycled marker properties.

# Arguments
- `model`: Model data container
- `outside_indices`: Indices of markers outside the model domain

"""
function reset_recycled_marker_properties!(
    model::ModelData,
    outside_indices::Vector{Int64}
)::Nothing
    nrecycle = length(outside_indices)
    matid_sticky_air = model.materials.dicts.matid_types["StickyAir"][1]
    temperature_reset = model.bcs.parameters.temperature.temperature_top.value

    for ioutside in 1:nrecycle
        imarker = outside_indices[ioutside]
        reset_marker!(model, imarker, matid_sticky_air, temperature_reset=temperature_reset)
    end
    return nothing
end

"""
    calculate_recycled_locations_along_top_boundary(
        nrecycle::Int,
        model::ModelData
    )::Tuple{Vector{Float64}, Vector{Float64}}

Calculate recycled locations.

# Returns
- `newx::Vector{Float64}`: New recycled x-locations.
- `newy::Vector{Float64}`: New recycled y-locations.
"""
function calculate_recycled_locations_along_top_boundary(
    nrecycle::Int,
    model::ModelData
)::Tuple{Vector{Float64}, Vector{Float64}}
    timestep = model.timestep.parameters.main_time_loop.timestep.value
    velocity_y_upper = model.bcs.parameters.velocity.vyu.value
    xsize = model.grids.parameters.geometry.xsize.value

    displacement = velocity_y_upper * timestep
    x_initial = 0.0
    newx = zeros(nrecycle)
    newy = zeros(nrecycle)

    for j in 1:nrecycle
        random_num = rand()
        newy[j] = random_num * displacement
        newx[j] = x_initial + random_num * xsize
    end
    return newx, newy
end

end # module