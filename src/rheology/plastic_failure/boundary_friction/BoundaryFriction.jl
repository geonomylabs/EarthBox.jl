module BoundaryFriction

import EarthBox.ModelDataContainer: ModelData
import EarthBox.GridFuncs: check_in_domain
export apply_boundary_friction, calculate_properties

function apply_boundary_friction(model::ModelData)::Nothing
    # Unpack data structures for better performance
    boundary_friction = model.materials.parameters.boundary_friction
    boundary_friction_width = boundary_friction.boundary_friction_width.value
    boundary_friction_angle = boundary_friction.boundary_friction_angle.value
    boundary_cohesion = boundary_friction.boundary_cohesion.value

    timesum = model.timestep.parameters.main_time_loop.timesum.value
    ysize = model.grids.parameters.geometry.ysize.value

    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array

    marker_cohesion = model.markers.arrays.rheology.marker_cohesion.array
    marker_fric = model.markers.arrays.rheology.marker_fric.array

    velocity = model.bcs.parameters.velocity
    velocity_internal_x = velocity.velocity_internal_x.value

    mobile_wall = model.geometry.parameters.mobile_wall
    y_top_mobile_wall = mobile_wall.y_top_mobile_wall.value
    x_left_mobile_wall = mobile_wall.x_left_mobile_wall.value

    marker_matid = model.markers.arrays.material.marker_matid.array
    matid_domains_dict = model.materials.dicts.matid_domains
    matid_sticky_air = matid_domains_dict["Atmosphere"]
    matid_pdms_layer = matid_domains_dict["PDMSLayer"]

    geometry = model.grids.parameters.geometry
    marknum = model.markers.parameters.distribution.marknum.value

    Threads.@threads for imarker in 1:marknum
        x_marker = marker_x[imarker]
        y_marker = marker_y[imarker]
        if check_in_domain(geometry, x_marker, y_marker)
            cohesion = marker_cohesion[imarker]
            sine_friction_angle = marker_fric[imarker]
            matid = marker_matid[imarker]
            
            (
                in_boundary_friction_zone, cohesion_new, sine_friction_angle_new
            ) = calculate_properties(
                cohesion, sine_friction_angle, timesum, ysize,
                x_left_mobile_wall, y_top_mobile_wall, boundary_friction_width,
                velocity_internal_x, x_marker, y_marker,
                boundary_friction_angle, boundary_cohesion, matid,
                matid_sticky_air, matid_pdms_layer
            )
            
            if in_boundary_friction_zone
                cohesion = cohesion_new
                sine_friction_angle = sine_friction_angle_new
            end
            marker_cohesion[imarker] = cohesion
            marker_fric[imarker] = sine_friction_angle
        end
    end
    return nothing
end

"""
Apply boundary friction cohesion and friction angle.

# Arguments
- `cohesion::Float64`: Marker cohesion (Pa)
- `sine_friction_angle::Float64`: Sine of friction angle (Degrees)
- `timesum::Float64`: Model time in seconds
- `ysize::Float64`: Model ysize in meters
- `x_left_mobile_wall::Float64`: Left edge of mobile wall (meters)
- `y_top_mobile_wall::Float64`: Top edge of mobile wall (meters)
- `boundary_friction_width::Float64`: Width of boundary friction zone (meters)
- `velocity_internal_x::Float64`: Prescribed shortening velocity (m/s)
- `x_marker::Float64`: Marker x-location (meters)
- `y_marker::Float64`: Marker y-location (meters)
- `boundary_friction_angle::Float64`: Friction angle of boundary friction zone (degrees)
- `boundary_cohesion::Float64`: Cohesion of boundary friction zone (Pa)
- `matid::Float64`: Material ID of marker
- `matid_sticky_air::Int`: Material ID of sticky air
- `matid_pdms_layer::Int`: Material ID of PDMS layer

"""
function calculate_properties(
    cohesion::Float64,
    sine_friction_angle::Float64,
    timesum::Float64,
    ysize::Float64,
    x_left_mobile_wall::Float64,
    y_top_mobile_wall::Float64,
    boundary_friction_width::Float64,
    velocity_internal_x::Float64,
    x_marker::Float64,
    y_marker::Float64,
    boundary_friction_angle::Float64,
    boundary_cohesion::Float64,
    matid::Int16,
    matid_sticky_air::Int16,
    matid_pdms_layer::Int16
)::Tuple{Bool, Float64, Float64}
    
    x_left_mobile_wall = x_left_mobile_wall + velocity_internal_x * timesum
    x_left_boundary_friction = x_left_mobile_wall - boundary_friction_width

    in_boundary_friction_zone = false
    
    # Vertical boundary friction zone along left edge of model
    if x_marker < boundary_friction_width
        in_boundary_friction_zone = true
    # Horizontal boundary friction zone along bottom edge of model
    elseif y_marker > ysize - boundary_friction_width
        in_boundary_friction_zone = true
    # Boundary friction zone along left edge of mobile wall
    elseif x_left_boundary_friction < x_marker < x_left_mobile_wall
        if y_marker > y_top_mobile_wall
            in_boundary_friction_zone = true
        end
    end
    
    # Don't use boundary friction in sticky air or PDMS layer
    if matid == matid_sticky_air || matid == matid_pdms_layer
        in_boundary_friction_zone = false
    end

    if in_boundary_friction_zone
        # Set cohesion to zero and friction angle to 19 degrees
        cohesion_new = boundary_cohesion
        sine_friction_angle_new = sin(boundary_friction_angle / 180.0 * π)
    else
        # Do not modify cohesion and friction angle
        cohesion_new = cohesion
        sine_friction_angle_new = sine_friction_angle
    end
    
    return (in_boundary_friction_zone, cohesion_new, sine_friction_angle_new)
end

end # module 