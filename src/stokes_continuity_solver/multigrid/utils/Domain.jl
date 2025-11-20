module Domain

function is_pressure_node(
    i::Int64, 
    j::Int64, 
    k::Int64, 
    ynum::Int64,
    xnum::Int64,
    znum::Int64
)::Bool
    return i < ynum && j < xnum && k < znum
end

function is_pressure_node_2d(
    i::Int64, 
    j::Int64, 
    ynum::Int64,
    xnum::Int64
)::Bool
    return i < ynum && j < xnum
end

function is_shearxy_node(
    k::Int64,
    znum::Int64
)::Bool
    return k < znum
end

function is_shearxz_node(
    i::Int64,
    ynum::Int64
)::Bool
    return i < ynum
end

function is_shearyz_node(
    j::Int64,
    xnum::Int64
)::Bool
    return j < xnum
end

function is_associated_with_internal_pressure_node(
    i::Int64, 
    j::Int64, 
    k::Int64, 
    ynum::Int64
)::Bool
    return i < ynum && j > 1 && k > 1
end

function is_associated_with_internal_pressure_node_2d(
    i::Int64, 
    j::Int64, 
    ynum::Int64
)::Bool
    return i < ynum && j > 1
end

function is_right_part_of_x_stokes_equation(
    i::Int64, 
    j::Int64, 
    k::Int64, 
    xnum::Int64
)::Bool
    return j > 1 && i > 1 && k > 1 && j < xnum
end

function is_right_part_of_x_stokes_equation_2d(
    i::Int64, 
    j::Int64, 
    xnum::Int64
)::Bool
    return 1 < j < xnum && i > 1
end

function is_right_part_of_y_stokes_equation(
    i::Int64, 
    j::Int64, 
    k::Int64, 
    ynum::Int64
)::Bool
    return j > 1 && k > 1 && i > 1 && i < ynum
end

function is_right_part_of_y_stokes_equation_2d(
    i::Int64, 
    j::Int64, 
    ynum::Int64
)::Bool
    return 1 < i < ynum && j > 1
end

function is_right_part_of_z_stokes_equation(
    i::Int64, 
    j::Int64, 
    k::Int64, 
    znum::Int64
)::Bool
    return j > 1 && i > 1 && k > 1 && k < znum
end

function on_vx_boundary3d(
    i::Int64, 
    j::Int64, 
    k::Int64,
    ynum::Int64,
    xnum::Int64,
    znum::Int64
)::Bool
    ynum_vx = get_ynum_vx_grid(ynum)
    znum_vx = get_znum_vx_grid(znum)
    return i == 1 || i == ynum_vx || j == 1 || j == xnum || k == 1 || k == znum_vx
end

function on_vx_boundary2d(
    i::Int64, 
    j::Int64, 
    ynum::Int64,
    xnum::Int64
)::Bool
    ynum_vx = get_ynum_vx_grid(ynum)
    return i == 1 || i == ynum_vx || j == 1 || j == xnum
end

function get_ynum_vx_grid(ynum::Int64)::Int64
    return ynum + 1
end

function get_znum_vx_grid(znum::Int64)::Int64
    return znum + 1
end

function on_vy_boundary3d(
    i::Int64, 
    j::Int64, 
    k::Int64,
    ynum::Int64,
    xnum::Int64,
    znum::Int64
)::Bool
    xnum_vy = get_xnum_vy_grid(xnum)
    znum_vy = get_znum_vy_grid(znum)
    return i == 1 || i == ynum || j == 1 || j == xnum_vy || k == 1 || k == znum_vy
end

function on_vy_boundary2d(
    i::Int64, 
    j::Int64, 
    ynum::Int64,
    xnum::Int64
)::Bool
    xnum_vy = get_xnum_vy_grid(xnum)
    return i == 1 || i == ynum || j == 1 || j == xnum_vy
end

function get_xnum_vy_grid(xnum::Int64)::Int64
    return xnum + 1
end

function get_znum_vy_grid(znum::Int64)::Int64
    return znum + 1
end

function on_vz_boundary3d(
    i::Int64, 
    j::Int64, 
    k::Int64,
    ynum::Int64,
    xnum::Int64,
    znum::Int64
)::Bool
    xnum_vz = get_xnum_vz_grid(xnum)
    ynum_vz = get_ynum_vz_grid(ynum)
    return i == 1 || i == ynum_vz || j == 1 || j == xnum_vz || k == 1 || k == znum
end

function get_xnum_vz_grid(xnum::Int64)::Int64
    return xnum + 1
end

function get_ynum_vz_grid(ynum::Int64)::Int64
    return ynum + 1
end

function get_ynum_pressure_grid(ynum::Int64)::Int64
    return ynum - 1
end

function get_xnum_pressure_grid(xnum::Int64)::Int64
    return xnum - 1
end

function get_znum_pressure_grid(znum::Int64)::Int64
    return znum - 1
end

end # module