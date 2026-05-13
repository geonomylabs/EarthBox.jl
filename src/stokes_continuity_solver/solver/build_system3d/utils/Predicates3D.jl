"""
Local predicates needed by the 3D Stokes-continuity build system that are
not yet available from `EarthBox.Domain`. These cover:

  - z-axis (front/back) versions of the basic-grid boundary checks.
  - 3D versions of the prescribed-internal-velocity-zone tests for
    x-Stokes, y-Stokes and z-Stokes.
  - 3D versions of the "found internal cell" tests for the four equations.

The x- and y-axis 2D predicates already exported by `EarthBox.Domain`
(`basic_cell_along_left_boundary`, `basic_cell_along_upper_boundary`,
`basic_cell_in_last_two_rightmost_columns`, etc.) are reused unchanged.

`bintern_zone` and `bintern_velocity` follow a 3D layout that extends the
2D convention. Each prescribed-velocity zone uses six consecutive entries:

  bintern_zone[1..6]  : vx zone — [j_col, min_i, max_i, min_k, max_k, unused]
                        (no prescription if bintern_zone[1] is negative)
  bintern_zone[7..12] : vy zone — [i_row, min_j, max_j, min_k, max_k, unused]
                        (no prescription if bintern_zone[7] is negative)
  bintern_zone[13..18]: vz zone — [k_slab, min_j, max_j, min_i, max_i, unused]
                        (no prescription if bintern_zone[13] is negative)

  bintern_velocity[1] : prescribed vx velocity (m/s)
  bintern_velocity[2] : prescribed vy velocity (m/s)
  bintern_velocity[3] : prescribed vz velocity (m/s)

Upstream is responsible for sizing these arrays once the 3D wiring is
implemented; this module only consumes them.
"""
module Predicates3D

export basic_cell_along_front_boundary
export basic_cell_along_back_boundary
export basic_cell_in_last_two_back_slabs
export basic_node_on_front_boundary
export basic_node_on_back_boundary
export found_internal_cell_x_stokes_3d
export found_internal_cell_y_stokes_3d
export found_internal_cell_z_stokes_3d
export found_internal_cell_continuity_3d
export outside_prescribed_velocity_zone_x_stokes_3d
export outside_prescribed_velocity_zone_y_stokes_3d
export outside_prescribed_velocity_zone_z_stokes_3d

@inline basic_cell_along_front_boundary(k::Int)::Bool = k == 1

@inline basic_cell_along_back_boundary(k::Int, znum::Int)::Bool = k == znum - 1

@inline basic_cell_in_last_two_back_slabs(k::Int, znum::Int)::Bool = k >= znum - 2

@inline basic_node_on_front_boundary(k::Int)::Bool = k == 1

@inline basic_node_on_back_boundary(k::Int, znum::Int)::Bool = k == znum

"""
    found_internal_cell_x_stokes_3d(j, xnum)

x-Stokes is applied at internal cells, which excludes cells in the
rightmost basic column (where vxC is prescribed).
"""
@inline found_internal_cell_x_stokes_3d(j::Int, xnum::Int)::Bool = j < xnum - 1

"""
    found_internal_cell_y_stokes_3d(i, ynum)

y-Stokes is applied at internal cells, which excludes cells in the
bottommost basic row (where vyC is prescribed).
"""
@inline found_internal_cell_y_stokes_3d(i::Int, ynum::Int)::Bool = i < ynum - 1

"""
    found_internal_cell_z_stokes_3d(k, znum)

z-Stokes is applied at internal cells, which excludes cells in the
backmost basic slab (where vzC is prescribed).
"""
@inline found_internal_cell_z_stokes_3d(k::Int, znum::Int)::Bool = k < znum - 1

"""
    found_internal_cell_continuity_3d(i, j, k, xnum, ynum, znum, mode)

A continuity cell is internal unless it carries a pressure boundary
condition. The mode parameter selects which cells are flagged as
boundary cells:

  * mode = 0: only the upper-left-front corner cell (i=j=k=1).
  * mode = 1: cells on the top (i=1) and bottom (i=ynum-1) y-faces.
  * mode = 2: cells on the left (j=1) and right (j=xnum-1) x-faces.
  * mode = 3: cells on the front (k=1) and back (k=znum-1) z-faces (new
              in 3D).
"""
@inline function found_internal_cell_continuity_3d(
    i::Int, j::Int, k::Int,
    xnum::Int, ynum::Int, znum::Int,
    mode::Int
)::Bool
    if mode == 0
        return !(i == 1 && j == 1 && k == 1)
    elseif mode == 1
        return i != 1 && i != ynum - 1
    elseif mode == 2
        return j != 1 && j != xnum - 1
    elseif mode == 3
        return k != 1 && k != znum - 1
    end
    return true
end

"""
    outside_prescribed_velocity_zone_x_stokes_3d(i, j, k, bintern_zone)

True iff cell (i, j, k) lies outside the prescribed-vx zone. The zone
is described by `bintern_zone[1..6]`; if `bintern_zone[1] < 0` no
prescription is active and every cell is considered "outside".
"""
@inline function outside_prescribed_velocity_zone_x_stokes_3d(
    i::Int, j::Int, k::Int, bintern_zone::Vector{Int64}
)::Bool
    bintern_zone[1] < 0 && return true
    inside = (
        j == bintern_zone[1] &&
        i >= bintern_zone[2] && i <= bintern_zone[3] &&
        k >= bintern_zone[4] && k <= bintern_zone[5]
    )
    return !inside
end

"""
    outside_prescribed_velocity_zone_y_stokes_3d(i, j, k, bintern_zone)

True iff cell (i, j, k) lies outside the prescribed-vy zone described by
`bintern_zone[7..12]`. No prescription is active if `bintern_zone[7] < 0`.
"""
@inline function outside_prescribed_velocity_zone_y_stokes_3d(
    i::Int, j::Int, k::Int, bintern_zone::Vector{Int64}
)::Bool
    bintern_zone[7] < 0 && return true
    inside = (
        i == bintern_zone[7] &&
        j >= bintern_zone[8] && j <= bintern_zone[9] &&
        k >= bintern_zone[10] && k <= bintern_zone[11]
    )
    return !inside
end

"""
    outside_prescribed_velocity_zone_z_stokes_3d(i, j, k, bintern_zone)

True iff cell (i, j, k) lies outside the prescribed-vz zone described by
`bintern_zone[13..18]`. No prescription is active if `bintern_zone[13] < 0`.
"""
@inline function outside_prescribed_velocity_zone_z_stokes_3d(
    i::Int, j::Int, k::Int, bintern_zone::Vector{Int64}
)::Bool
    bintern_zone[13] < 0 && return true
    inside = (
        k == bintern_zone[13] &&
        j >= bintern_zone[14] && j <= bintern_zone[15] &&
        i >= bintern_zone[16] && i <= bintern_zone[17]
    )
    return !inside
end

end # module
