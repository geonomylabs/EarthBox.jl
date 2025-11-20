module SerpentineRheology

import EarthBox.ModelDataContainer: ModelData
import EarthBox: GridFuncs

""" Update marker flow viscosity for serpentinization.

# Arguments
- `model::ModelData`: Main model interface object
- `mantle_ids::Vector{Int16}`: Mantle material ids
- `max_strength_reduction_factor::Float64`: Maximum strength reduction factor

# Updated Arrays
- `model.markers.arrays.rheology.marker_eta_flow`: Marker effective flow viscosity (Pa.s)
"""
function update_marker_flow_viscosity!(
    model::ModelData,
    mantle_ids::Vector{Int16},
    max_strength_reduction_factor::Float64=10.0
)::Nothing
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marker_eta_flow = model.markers.arrays.rheology.marker_eta_flow.array
    marker_matid = model.markers.arrays.material.marker_matid.array
    marker_serpentinization = model.markers.arrays.material.marker_serpentinization.array
    viscosity_minimum = model.materials.parameters.viscosity_limits.viscosity_min
    viscosity_maximum = model.materials.parameters.viscosity_limits.viscosity_max
    geometry = model.grids.parameters.geometry
    marknum = model.markers.parameters.distribution.marknum

    Threads.@threads for imarker in 1:marknum
        x_marker = marker_x[imarker]
        y_marker = marker_y[imarker]
        if GridFuncs.check_in_domain(geometry, x_marker, y_marker)
            serpentinization_ratio = marker_serpentinization[imarker]
            mat_id = marker_matid[imarker]
            if mat_id in mantle_ids && serpentinization_ratio > 0.0
                viscosity_flow = marker_eta_flow[imarker]
                reduction_factor = calculate_strength_reduction_factor(
                    serpentinization_ratio,
                    max_strength_reduction_factor
                )
                viscosity_flow = viscosity_flow / reduction_factor
                viscosity_flow = apply_viscosity_limits(
                    viscosity_minimum,
                    viscosity_maximum,
                    viscosity_flow
                )
                marker_eta_flow[imarker] = viscosity_flow
            end
        end
    end
    return nothing
end

""" Limit effective viscosity for numerical stability.
"""
function apply_viscosity_limits(
    viscosity_minimum::Float64,
    viscosity_maximum::Float64,
    viscosity::Float64
)::Float64
    if viscosity < viscosity_minimum
        viscosity = viscosity_minimum
    elseif viscosity > viscosity_maximum
        viscosity = viscosity_maximum
    end
    return viscosity
end

""" Update marker cohesion and friction coefficient arrays for serpentinization.

# Arguments
- `model::ModelData`: Main model interface object
- `mantle_ids::Vector{Int16}`: Mantle material ids
- `max_strength_reduction_factor::Float64`: Maximum strength reduction factor

# Updated Arrays
- `model.markers.arrays.rheology.marker_cohesion`: Marker cohesion (Pa)
- `model.markers.arrays.rheology.marker_fric`: Marker sine of friction angle
"""
function update_marker_failure_properties!(
    model::ModelData,
    mantle_ids::Vector{Int16},
    max_strength_reduction_factor::Float64=2.0
)
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marker_matid = model.markers.arrays.material.marker_matid.array
    marker_cohesion = model.markers.arrays.rheology.marker_cohesion.array
    marker_fric = model.markers.arrays.rheology.marker_fric.array
    marker_serpentinization = model.markers.arrays.material.marker_serpentinization.array
    geometry = model.grids.parameters.geometry
    marknum = model.markers.parameters.distribution.marknum

    Threads.@threads for imarker in 1:marknum
        x_marker = marker_x[imarker]
        y_marker = marker_y[imarker]
        if GridFuncs.check_in_domain(geometry, x_marker, y_marker)
            serpentinization_ratio = marker_serpentinization[imarker]
            mat_id = marker_matid[imarker]
            if mat_id in mantle_ids && serpentinization_ratio > 0.0
                cohesion = marker_cohesion[imarker]
                sin_friction_angle = marker_fric[imarker]
                reduction_factor = calculate_strength_reduction_factor(
                    serpentinization_ratio,
                    max_strength_reduction_factor
                )
                marker_cohesion[imarker] = cohesion / reduction_factor
                marker_fric[imarker] = sin_friction_angle / reduction_factor
            end
        end
    end
end

""" Calculate strength reduction factor for serpentinization.

# Arguments
- `serpentinization_ratio::Float64`: Serpentinization ratio
- `max_strength_reduction_factor::Float64`: Maximum strength reduction factor

# Returns
- `reduction_factor::Float64`: Strength reduction factor
"""
function calculate_strength_reduction_factor(
    serpentinization_ratio::Float64,
    max_strength_reduction_factor::Float64
)::Float64
    return 1.0 + (max_strength_reduction_factor - 1.0) * serpentinization_ratio
end

end