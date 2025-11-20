module MobileWall

using EarthBox.ModelDataContainer: ModelData

"""
    in_mobile_wall(
        model::ModelData, # type: ignore
        i::Int,
        j::Int
    ) -> Tuple{Bool, Float64, Float64}

Check if marker is in mobile wall.

# Arguments
- `model`: Model data structure
- `i`: Row index of marker
- `j`: Column index of marker

# Returns
- Tuple containing (is_in_wall::Bool, x_node_m::Float64, y_node_m::Float64)
"""
function in_mobile_wall(
    model::ModelData, i::Int, j::Int
)::Tuple{Bool, Float64, Float64}
    gridx = model.grids.arrays.basic.gridx_b.array
    gridy = model.grids.arrays.basic.gridy_b.array
    xsize = model.grids.parameters.geometry.xsize.value
    ysize = model.grids.parameters.geometry.ysize.value
    
    x_left_mobile_wall = xsize - 0.01
    x_right_mobile_wall = xsize - 0.005
    y_top_mobile_wall = 0.005
    y_bottom_mobile_wall = ysize - 0.002
    
    x_node_m = gridx[j]
    y_node_m = gridy[i]
    
    check_flag = define_mobile_wall(
        0, 1, x_node_m, y_node_m,
        x_left_mobile_wall, x_right_mobile_wall,
        y_top_mobile_wall, y_bottom_mobile_wall
    )
    
    check = check_flag == 1
    return check, x_node_m, y_node_m
end

"""
    in_mobile_wall_xy(model::ModelData, xm::Float64, ym::Float64) -> Bool

Check if marker is in mobile wall using x,y coordinates.

# Arguments
- `model`: Model data structure
- `xm`: X-coordinate of marker
- `ym`: Y-coordinate of marker

# Returns
- Bool indicating if marker is in mobile wall
"""
function in_mobile_wall_xy(model::ModelData, xm::Float64, ym::Float64)
    xsize = model.grids.parameters.geometry.xsize.value
    ysize = model.grids.parameters.geometry.ysize.value
    
    x_left_mobile_wall = xsize - 0.01
    x_right_mobile_wall = xsize - 0.005
    y_top_mobile_wall = 0.005
    y_bottom_mobile_wall = ysize - 0.002
    
    check_flag = define_mobile_wall(
        0, 1, xm, ym,
        x_left_mobile_wall, x_right_mobile_wall,
        y_top_mobile_wall, y_bottom_mobile_wall
    )
    
    return check_flag == 1
end

"""
    in_mobile_wall_transient(
        x_node::Float64,
        y_node::Float64,
        x_left_mobile_wall::Float64,
        x_right_mobile_wall::Float64,
        y_top_mobile_wall::Float64,
        y_bottom_mobile_wall::Float64, 
        velocity_internal_x::Float64,
        timesum::Float64
    ) -> Bool

Check if marker is in the mobile wall with transient movement.

# Arguments
- `x_node`: X-coordinate of marker (m)
- `y_node`: Y-coordinate of marker (m)
- `x_left_mobile_wall`: X-coordinate of left mobile wall (m)
- `x_right_mobile_wall`: X-coordinate of right mobile wall (m)
- `y_top_mobile_wall`: Y-coordinate of top mobile wall (m)
- `y_bottom_mobile_wall`: Y-coordinate of bottom mobile wall (m)
- `velocity_internal_x`: X-velocity of mobile wall (m/s)
- `timesum`: Model time in seconds

# Returns
- Bool indicating if marker is in mobile wall
"""
@inline function in_mobile_wall_transient(
    x_node::Float64, 
    y_node::Float64, 
    x_left_mobile_wall::Float64,
    x_right_mobile_wall::Float64, 
    y_top_mobile_wall::Float64,
    y_bottom_mobile_wall::Float64, 
    velocity_internal_x::Float64,
    timesum::Float64
)::Bool
    x_left_mobile_wall += velocity_internal_x * timesum
    x_right_mobile_wall += velocity_internal_x * timesum
    return is_in_mobile_wall(
        x_node, y_node, x_left_mobile_wall, x_right_mobile_wall,
        y_top_mobile_wall, y_bottom_mobile_wall
    )
end

"""
    define_mobile_wall(
        matid::Int16, 
        matid_mobile_wall::Int16, 
        x_marker::Float64,
        y_marker::Float64, 
        x_left_mobile_wall::Float64,
        x_right_mobile_wall::Float64,
        y_top_mobile_wall::Float64,
        y_bottom_mobile_wall::Float64
    ) -> Int16

Define mobile wall material ID.

# Arguments
- `matid`: Material ID of marker
- `matid_mobile_wall`: Material ID of mobile wall
- `x_marker`: X-coordinate of marker
- `y_marker`: Y-coordinate of marker
- `x_left_mobile_wall`: X-coordinate of left mobile wall
- `x_right_mobile_wall`: X-coordinate of right mobile wall
- `y_top_mobile_wall`: Y-coordinate of top mobile wall
- `y_bottom_mobile_wall`: Y-coordinate of bottom mobile wall

# Returns
- Material ID of marker
"""
function define_mobile_wall(
    matid::Int16, 
    matid_mobile_wall::Int16, 
    x_marker::Float64,
    y_marker::Float64,
    x_left_mobile_wall::Float64,
    x_right_mobile_wall::Float64, 
    y_top_mobile_wall::Float64,
    y_bottom_mobile_wall::Float64
)::Int16
    if is_in_mobile_wall(
        x_marker, y_marker, x_left_mobile_wall, x_right_mobile_wall,
        y_top_mobile_wall, y_bottom_mobile_wall
    )
        return matid_mobile_wall
    end
    return matid
end

"""
    is_in_mobile_wall(
        x_marker::Float64, 
        y_marker::Float64, 
        x_left_mobile_wall::Float64,
        x_right_mobile_wall::Float64, 
        y_top_mobile_wall::Float64,
        y_bottom_mobile_wall::Float64
    ) -> Bool

Check if marker is in mobile wall.

# Arguments
- `x_marker`: X-coordinate of marker
- `y_marker`: Y-coordinate of marker
- `x_left_mobile_wall`: X-coordinate of left mobile wall
- `x_right_mobile_wall`: X-coordinate of right mobile wall
- `y_top_mobile_wall`: Y-coordinate of top mobile wall
- `y_bottom_mobile_wall`: Y-coordinate of bottom mobile wall

# Returns
- Bool indicating if marker is in mobile wall
"""
function is_in_mobile_wall(
    x_marker::Float64, 
    y_marker::Float64, 
    x_left_mobile_wall::Float64,
    x_right_mobile_wall::Float64, 
    y_top_mobile_wall::Float64,
    y_bottom_mobile_wall::Float64
)::Bool
    return (x_left_mobile_wall <= x_marker <= x_right_mobile_wall) &&
           (y_top_mobile_wall <= y_marker <= y_bottom_mobile_wall)
end

"""
    define_plate_extension(
        matid::Int16, 
        x_marker::Float64, 
        y_marker::Float64,
        matid_plate_extension::Int16, 
        plate_extension_thickness::Float64,
        plate_extension_width::Float64, 
        x_left_mobile_wall::Float64,
        y_bottom_mobile_wall::Float64
    ) -> Int16

Define plate extension material ID.

# Arguments
- `matid`: Material ID of marker
- `x_marker`: X-coordinate of marker
- `y_marker`: Y-coordinate of marker
- `matid_plate_extension`: Material ID of plate extension
- `plate_extension_thickness`: Thickness of plate extension
- `plate_extension_width`: Width of plate extension
- `x_left_mobile_wall`: X-coordinate of left mobile wall
- `y_bottom_mobile_wall`: Y-coordinate of bottom mobile wall

# Returns
- Material ID of marker
"""
function define_plate_extension(
    matid::Int16, 
    x_marker::Float64, 
    y_marker::Float64,
    matid_plate_extension::Int16,
    plate_extension_thickness::Float64,
    plate_extension_width::Float64, 
    x_left_mobile_wall::Float64,
    y_bottom_mobile_wall::Float64
)::Int16
    xmin = x_left_mobile_wall - plate_extension_width
    xmax = x_left_mobile_wall
    ymin = y_bottom_mobile_wall - plate_extension_thickness
    ymax = y_bottom_mobile_wall
    if is_in_plate_extension(x_marker, y_marker, xmin, xmax, ymin, ymax)
        return matid_plate_extension
    end
    return matid
end

"""
    in_plate_extension_transient(
        x_node::Float64,
        y_node::Float64,
        plate_extension_thickness::Float64,
        plate_extension_width::Float64,
        x_left_mobile_wall::Float64,
        y_bottom_mobile_wall::Float64,
        velocity_internal_x::Float64,
        timesum::Float64
    ) -> Bool

Check if marker is in the plate extension with transient movement.

# Arguments
- `x_node`: X-coordinate of marker
- `y_node`: Y-coordinate of marker
- `plate_extension_thickness`: Thickness of plate extension
- `plate_extension_width`: Width of plate extension
- `x_left_mobile_wall`: X-coordinate of left mobile wall
- `y_bottom_mobile_wall`: Y-coordinate of bottom mobile wall
- `velocity_internal_x`: X-velocity of mobile wall
- `timesum`: Model time in seconds

# Returns
- Bool indicating if marker is in plate extension
"""
@inline function in_plate_extension_transient(
    x_node::Float64, 
    y_node::Float64,
    plate_extension_thickness::Float64,
    plate_extension_width::Float64,
    x_left_mobile_wall::Float64,
    y_bottom_mobile_wall::Float64,
    velocity_internal_x::Float64,
    timesum::Float64
)::Bool
    x_left_mobile_wall += velocity_internal_x * timesum
    xmin = x_left_mobile_wall - plate_extension_width
    xmax = x_left_mobile_wall
    ymin = y_bottom_mobile_wall - plate_extension_thickness
    ymax = y_bottom_mobile_wall
    return is_in_plate_extension(x_node, y_node, xmin, xmax, ymin, ymax)
end

"""
    is_in_plate_extension(
        x_node::Float64, 
        y_node::Float64, 
        xmin::Float64,
        xmax::Float64, 
        ymin::Float64, 
        ymax::Float64
    ) -> Bool

Check if marker is in the plate extension.

# Arguments
- `x_node`: X-coordinate of marker
- `y_node`: Y-coordinate of marker
- `xmin`: Minimum x-coordinate of plate extension
- `xmax`: Maximum x-coordinate of plate extension
- `ymin`: Minimum y-coordinate of plate extension
- `ymax`: Maximum y-coordinate of plate extension

# Returns
- Bool indicating if marker is in plate extension
"""
@inline function is_in_plate_extension(
    x_node::Float64,
    y_node::Float64,
    xmin::Float64,
    xmax::Float64,
    ymin::Float64,
    ymax::Float64
)::Bool
    return (xmin <= x_node <= xmax) && (ymin <= y_node <= ymax)
end

end # module MobileWall 