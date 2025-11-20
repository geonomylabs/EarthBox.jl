module Domain

""" This was referred to as on_boundary() in old heat domain module
"""
@inline function basic_node_on_boundary(i::Int, j::Int, xnum::Int, ynum::Int)::Bool
    return j == 1 || j == xnum || i == 1 || i == ynum
end

""" This was referred to as on_left_boundary() in old heat domain module
"""
function basic_node_on_left_boundary(j::Int)::Bool
    return j == 1
end

""" This was referred to as on_right_boundary() in old heat domain module
"""
function basic_node_on_right_boundary(j::Int, xnum::Int)::Bool
    return j == xnum
end

""" This was referred to as on_upper_boundary() in old heat domain module
"""
function basic_node_on_upper_inside_boundary(i::Int, j::Int, xnum::Int)::Bool
    return i == 1 && 1 < j < xnum
end

""" This was referred to as on_lower_boundary() in old heat domain module
"""
function basic_node_on_lower_inside_boundary(i::Int, j::Int, xnum::Int, ynum::Int)::Bool
    return i == ynum && 1 < j < xnum
end

function basic_node_along_lower_boundary(i::Int, ynum::Int)::Bool
    return i == ynum
end

function basic_cell_along_upper_boundary(i::Int)::Bool
    return i == 1
end

function basic_cell_along_lower_boundary(i::Int, ynum::Int)::Bool
    return i ≥ ynum - 1
end

function basic_cell_along_left_boundary(j::Int)::Bool
    return j == 1
end

function basic_cell_along_right_boundary(j::Int, xnum::Int)::Bool
    return j ≥ xnum - 1
end

function basic_cell_in_last_two_rightmost_columns(j::Int, xnum::Int)::Bool
    return j ≥ xnum - 2
end

function basic_cell_in_last_two_bottom_rows(i::Int, ynum::Int)::Bool
    return i ≥ ynum - 2
end

function outside_prescribed_velocity_zone_x_stokes(
    i::Int,
    j::Int,
    bintern_zone::Vector{Int64}
)::Bool
    check = false
    if (j != bintern_zone[1] || i < bintern_zone[2] || i > bintern_zone[3])
        check = true
    end
    return check
end

function outside_prescribed_velocity_zone_y_stokes(
    i::Int,
    j::Int,
    bintern_zone::Vector{Int64}
)::Bool
    check = false
    if (j != bintern_zone[5] || i < bintern_zone[6] || i > bintern_zone[7])
        check = true
    end
    return check
end

"""
    found_internal_cell_x_stokes(j::Int, xnum::Int)::Bool

Check if j of basic grid cell is two node steps from right boundary.

The first part of the conditional x-Stokes statement applies only to basic
grid cells with upper left nodes located at least two steps from the right
boundary of the basic grid. Consider the following visual:

  j=xnum-2       j=xnum-1      j=xnum
      o-------------x------------RB

where o refers to a node for which the first part of the conditional state
is applicable, x is a basic node and RB is a basic node located at at right
boundary of the basic grid.
"""
function found_internal_cell_x_stokes(j::Int, xnum::Int)::Bool
    check = false
    if j < xnum - 1
        check = true
    end
    return check
end

"""
    found_internal_cell_y_stokes(i::Int, ynum::Int)::Bool

Check if i of basic grid cell is two node steps from lower boundary.

The first part of this conditional statement applies only to basic cells
with upper left node located at least one steps from the bottom boundary
of the basic grid. Consider the following visual:

      o i=ynum-2
      |
      |
      |
      x i=ynum-1
      |
      |
      |
     LB i=ynum

where o refers to a node for which the first part of the conditional
statement is applicable, x refers to a non-applicable basic node and
LB is a basic node located at lower boundary of the basic grid.
"""
function found_internal_cell_y_stokes(i::Int, ynum::Int)::Bool
    check = false
    if i < ynum - 1
        check = true
    end
    return check
end

"""
    found_internal_cell_continuity(i::Int, j::Int, xnum::Int, ynum::Int, 
                                 pressure_bc_mode::Int)::Bool

Check if i, j of basic grid cell is an internal pressure cell.

Internal cells depend on pressure boundary condition type (pfirst_mode).
If pfirst_mode = 0, then the upper left corner cell is the boundary cell
while all others are internal cells for pressure. If pfirst_mode = 1 then
the top and bottom rows of cells are boundary cells while all others are
internal. If pfirst_mode = 2 the side columns of cells are boundary cells
while all others are internal.
"""
function found_internal_cell_continuity(
    i::Int,
    j::Int,
    xnum::Int,
    ynum::Int,
    pressure_bc_mode::Int
)::Bool
    check = false
    if pressure_bc_mode == 0
        if j > 1 || i > 1
            check = true
        end
    elseif pressure_bc_mode == 1
        if 1 < i < ynum - 1
            check = true
        end
    elseif pressure_bc_mode == 2
        if 1 < j < xnum - 1
            check = true
        end
    end
    return check
end

end # module 