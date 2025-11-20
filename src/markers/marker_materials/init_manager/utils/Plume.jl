module Plume

using ..InitManager.InitStructs: PlumeGeometry

"""
    in_plume_region_check(x_marker::Float64, y_marker::Float64,
                         plume_geometry::PlumeGeometry) -> Int

Check if marker is within plume head.

ihead_type = 0 = circular plume head
ihead_type = 1 = square plume head

# Arguments
- `x_marker`: X-coordinate of marker
- `y_marker`: Y-coordinate of marker
- `plume_geometry`: Plume geometry parameters

# Returns
- Int indicating if marker is in plume head (1 if true, 0 if false)
"""
function in_plume_region_check(
    x_marker::Float64, 
    y_marker::Float64,
    plume_geometry::PlumeGeometry
)::Int
    ihead_type = 0
    in_plume = false
    
    if ihead_type == 0
        # Circle
        dist = calculate_distance_from_plume_center(x_marker, y_marker, plume_geometry)
        if dist < plume_geometry.plume_radius
            in_plume = true
        end
    else
        # Square
        ymin, ymax, xmin, xmax = calculate_limits_of_square_plume_head(plume_geometry)
        if xmin < x_marker < xmax && ymin < y_marker < ymax
            in_plume = true
        end
    end
    
    return in_plume ? 1 : 0
end

"""
    calculate_distance_from_plume_center(x_marker::Float64, y_marker::Float64,
                                       plume_geometry::PlumeGeometry) -> Float64

Calculate distance from marker to plume head center.

# Arguments
- `x_marker`: X-coordinate of marker
- `y_marker`: Y-coordinate of marker
- `plume_geometry`: Plume geometry parameters

# Returns
- Distance from marker to plume-head center
"""
function calculate_distance_from_plume_center(
    x_marker::Float64, 
    y_marker::Float64,
    plume_geometry::PlumeGeometry
)::Float64
    plume_center_x = plume_geometry.plume_center_x
    plume_center_y = plume_geometry.plume_center_y
    dx_dist = x_marker - plume_center_x
    dy_dist = y_marker - plume_center_y
    return sqrt(dx_dist * dx_dist + dy_dist * dy_dist)
end

"""
    calculate_limits_of_square_plume_head(plume_geometry::PlumeGeometry) -> 
    Tuple{Float64, Float64, Float64, Float64}

Calculate limits of square plume head.

# Arguments
- `plume_geometry`: Plume geometry parameters

# Returns
- Tuple containing (ymin, ymax, xmin, xmax) limits
"""
function calculate_limits_of_square_plume_head(plume_geometry::PlumeGeometry)
    plume_radius = plume_geometry.plume_radius
    plume_center_x = plume_geometry.plume_center_x
    plume_center_y = plume_geometry.plume_center_y
    plume_head_thick = plume_geometry.plume_head_thick
    
    ymin = plume_center_y - plume_head_thick / 2.0
    ymax = plume_center_y + plume_head_thick / 2.0
    xmin = plume_center_x - plume_radius
    xmax = plume_center_x + plume_radius
    
    return ymin, ymax, xmin, xmax
end

end # module Plume 