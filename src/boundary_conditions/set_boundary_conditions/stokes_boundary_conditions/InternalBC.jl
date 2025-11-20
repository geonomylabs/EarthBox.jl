"""
    Boundary condition coefficients for internal prescribed velocity.
"""
module InternalBC

using EarthBox.ModelDataContainer: ModelData

"""
    turn_off_internal_velocity!(model::ModelData)::Nothing

Internal velocity option not used.

# Updated Arrays
- `model.bcs.arrays.internal.bintern.array::Vector{Int64}`: Array((8))
    - (1) Horizontal position of vx nodes with prescribed velocity [no such condition if negative]
    - (2) Min vertical position
    - (3) Max vertical position
    - (4) (not used)
    - (5) Horizontal position of vy nodes with prescribed velocity [no such condition if negative]
    - (6) Min vertical position
    - (7) Max vertical position
    - (8) (not used)

- `model.bcs.arrays.internal.bintern_v.array::Vector{Float64}`: Array((8))
    - (1) (not used)
    - (2) (not used)
    - (3) (not used)
    - (4) Prescribed shortening velocity, m/s
    - (5) (not used)
    - (6) (not used)
    - (7) (not used)
    - (8) Prescribed vertical velocity, m/s
"""
function turn_off_internal_velocity!(model::ModelData)::Nothing
    bintern_zone = model.bcs.arrays.internal.bintern_zone.array
    bintern_velocity = model.bcs.arrays.internal.bintern_velocity.array
    # If negative internal x-velocity bc's are ignored.
    bintern_zone[1] = -1
    bintern_zone[2] = 0
    bintern_zone[3] = 0
    bintern_zone[4] = 0.0
    bintern_velocity[4] = 0.0
    # if negative internal y-velocity bc's are ignored.
    bintern_zone[5] = -1
    bintern_zone[6] = 0
    bintern_zone[7] = 0
    bintern_zone[8] = 0.0
    bintern_velocity[8] = 0.0
    return nothing
end

"""
    set_internal_velocity_for_sandbox_with_mobile_wall(model::ModelData)

Prescribed internal velocity of "mobile wall" for sandbox experiments.

# Updated Arrays
- `model.bcs.arrays.internal.bintern_zone.array::Vector{Int64}`: Array((8))
    - (1) Horizontal position of vx nodes with prescribed velocity [no such condition if negative]
    - (2) Min vertical position
    - (3) Max vertical position
    - (4) (not used)
    - (5) Horizontal position of vy nodes with prescribed velocity [no such condition if negative]
    - (6) Min vertical position
    - (7) Max vertical position
    - (8) (not used)

- `model.bcs.arrays.internal.bintern_velocity.array::Vector{Float64}`: Array((8))
    - (1) (not used)
    - (2) (not used)
    - (3) (not used)
    - (4) Prescribed shortening velocity, m/s
    - (5) (not used)
    - (6) (not used)
    - (7) (not used)
    - (8) Prescribed vertical velocity, m/s
"""
function set_internal_velocity_for_sandbox_with_mobile_wall!(model::ModelData)::Nothing
    # Unpack model data
    bintern_zone = model.bcs.arrays.internal.bintern_zone.array
    bintern_velocity = model.bcs.arrays.internal.bintern_velocity.array

    velocity = model.bcs.parameters.velocity
    velocity_internal_x = velocity.velocity_internal_x.value
    velocity_internal_y = velocity.velocity_internal_y.value

    update_internal_velocity_zone_indices!(model)

    internal_zone = model.geometry.parameters.internal_velocity_zone

    xindex_vx_internal = internal_zone.xindex_vx_internal.value
    yindex_min_vx_internal = internal_zone.yindex_min_vx_internal.value
    yindex_max_vx_internal = internal_zone.yindex_max_vx_internal.value

    xindex_vy_internal = internal_zone.xindex_vy_internal.value
    yindex_min_vy_internal = internal_zone.yindex_min_vy_internal.value
    yindex_max_vy_internal = internal_zone.yindex_max_vy_internal.value

    #*************************
    # Internal x-velocity zone
    #*************************
    # Horizontal position of vx nodes with prescribed velocity.
    # If negative internal bc's are ignored.
    bintern_zone[1] = xindex_vx_internal # xnum - 9
    # Min vertical position
    bintern_zone[2] = yindex_min_vx_internal # 5
    # Max vertical position
    bintern_zone[3] = yindex_max_vx_internal # ynum - 4
    # Prescribed shortening velocity, m/s
    # bintern_zone[4] = velocity_internal_x # -2.5/100.0/3600.0 # Not used
    bintern_velocity[4] = velocity_internal_x # -2.5/100.0/3600.0
    #*************************
    # Internal y-velocity zone
    #*************************
    # Horizontal position of vy nodes with prescribed velocity.
    # If negative internal bc's are ignored.
    bintern_zone[5] = xindex_vy_internal # xnum - 9
    # Min vertical position
    bintern_zone[6] = yindex_min_vy_internal # ynum - 4
    # Max vertical position
    bintern_zone[7] = yindex_max_vy_internal # ynum - 4
    # Prescribed vertical velocity, m/s
    bintern_zone[8] = velocity_internal_y # 0.0
    bintern_velocity[8] = velocity_internal_y # 0.0
    return nothing
end

"""
    update_internal_velocity_zone_indices(model::ModelData)

Update internal velocity zone indices.

# Updated Parameters
- `model.geometry.parameters.internal_velocity_zone`
    - `xindex_vx_internal::Int64`: X-index of internal velocity x zone.
    - `yindex_min_vx_internal::Int64`: Minimum y-index of internal velocity x zone.
    - `yindex_max_vx_internal::Int64`: Maximum y-index of internal velocity x zone.
    - `xindex_vy_internal::Int64`: X-index of internal velocity y zone.
    - `yindex_min_vy_internal::Int64`: Minimum y-index of internal velocity y zone.
    - `yindex_max_vy_internal::Int64`: Maximum y-index of internal velocity y zone.
"""
function update_internal_velocity_zone_indices!(model::ModelData)::Nothing
    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value

    xsize = model.grids.parameters.geometry.xsize.value
    ysize = model.grids.parameters.geometry.ysize.value

    dx = xsize/float(xnum - 1)
    dy = ysize/float(ynum - 1)

    internal_zone = model.geometry.parameters.internal_velocity_zone

    x_vx_internal = internal_zone.x_vx_internal.value
    y_min_vx_internal = internal_zone.y_min_vx_internal.value
    y_max_vx_internal = internal_zone.y_max_vx_internal.value

    x_vy_internal = internal_zone.x_vy_internal.value
    y_min_vy_internal = internal_zone.y_min_vy_internal.value
    y_max_vy_internal = internal_zone.y_max_vy_internal.value

    xindex_vx_internal = Int(round(x_vx_internal/dx))
    yindex_min_vx_internal = Int(round(y_min_vx_internal/dy))
    yindex_max_vx_internal = Int(round(y_max_vx_internal/dy))

    internal_zone.xindex_vx_internal.value = xindex_vx_internal
    internal_zone.yindex_min_vx_internal.value = yindex_min_vx_internal
    internal_zone.yindex_max_vx_internal.value = yindex_max_vx_internal

    xindex_vy_internal = Int(round(x_vy_internal/dx))
    yindex_min_vy_internal = Int(round(y_min_vy_internal/dy))
    yindex_max_vy_internal = Int(round(y_max_vy_internal/dy))

    internal_zone.xindex_vy_internal.value = xindex_vy_internal
    internal_zone.yindex_min_vy_internal.value = yindex_min_vy_internal
    internal_zone.yindex_max_vy_internal.value = yindex_max_vy_internal
    return nothing
end

end # module InternalBC 