module VxBoundaryConditions

import ...LevelManager: LevelData, LevelData2d
import ...Domain: get_ynum_vx_grid, get_znum_vx_grid

function set_vx_boundary_conditions3d!(
    i::Int64,
    j::Int64,
    k::Int64,
    dims::NamedTuple,
    level_data::LevelData;
    bc_type::Symbol = :free_slip,
    vx_bc_left::Float64 = 0.0,
    vx_bc_right::Float64 = 0.0,
    vx_bc_upper::Float64 = 0.0,
    vx_bc_lower::Float64 = 0.0,
    vx_bc_back::Float64 = 0.0,
    vx_bc_front::Float64 = 0.0
)::Nothing
    vx = level_data.vx.array
    xnum = dims.xnum
    ynum_vx = dims.ynum_vx
    znum_vx = dims.znum_vx

    # Boundary condition for vx
    vx[i,j,k] = 0.0
    # Left boundary
    if j == 1
        if bc_type == :constant_velocity
            vx[i,j,k] = vx_bc_left
        end
    end
    # Right boundary
    if j == xnum
        if bc_type == :constant_velocity
            vx[i,j,k] = vx_bc_right
        end
    end
    # Upper Boundary
    if i == 1
        # Free slip
        if bc_type == :free_slip
            vx[i,j,k] = vx[i+1,j,k]
        # No slip
        elseif bc_type == :no_slip
            vx[i,j,k] = -vx[i+1,j,k]
        elseif bc_type == :constant_velocity
            vx[i,j,k] = vx_bc_upper
        end
    end
    # Lower boundary
    if i == ynum_vx
        # Free slip
        if bc_type == :free_slip
            vx[i,j,k] = vx[i-1,j,k]
        # No slip
        elseif bc_type == :no_slip
            vx[i,j,k] = -vx[i-1,j,k]
        elseif bc_type == :constant_velocity
            vx[i,j,k] = vx_bc_lower
        end
    end
    # Back Boundary
    if k == 1
        # Free slip
        if bc_type == :free_slip
            vx[i,j,k] = vx[i,j,k+1]
        # No slip
        elseif bc_type == :no_slip
            vx[i,j,k] = -vx[i,j,k+1]
        elseif bc_type == :constant_velocity
            vx[i,j,k] = vx_bc_back
        end
    end
    # Front boundary
    if k == znum_vx
        # Free slip
        if bc_type == :free_slip
            vx[i,j,k] = vx[i,j,k-1]
        # No slip
        elseif bc_type == :no_slip
            vx[i,j,k] = -vx[i,j,k-1]
        elseif bc_type == :constant_velocity
            vx[i,j,k] = vx_bc_front
        end
    end
    return nothing
end

function set_vx_boundary_conditions2d!(
    i::Int64,
    j::Int64,
    dims::NamedTuple,
    level_data::LevelData2d;
    bc_type::Symbol = :free_slip,
    vx_bc_left::Float64 = 0.0,
    vx_bc_right::Float64 = 0.0,
    vx_bc_upper::Float64 = 0.0,
    vx_bc_lower::Float64 = 0.0
)::Nothing
    vx = level_data.vx.array
    xnum = dims.xnum
    ynum_vx = dims.ynum_vx

    # Boundary condition for vx
    vx[i,j] = 0.0
    # Left boundary
    if j == 1
        if bc_type == :constant_velocity
            vx[i,j] = vx_bc_left
        end
    end
    # Right boundary
    if j == xnum
        if bc_type == :constant_velocity
            vx[i,j,k] = vx_bc_right
        end
    end
    # Upper Boundary
    if i == 1
        # Free slip
        if bc_type == :free_slip
            vx[i,j] = vx[i+1,j]
        # No slip
        elseif bc_type == :no_slip
            vx[i,j] = -vx[i+1,j]
        elseif bc_type == :constant_velocity
            vx[i,j] = vx_bc_upper
        end
    end
    # Lower boundary
    if i == ynum_vx
        # Free slip
        if bc_type == :free_slip
            vx[i,j] = vx[i-1,j]
        # No slip
        elseif bc_type == :no_slip
            vx[i,j] = -vx[i-1,j]
        elseif bc_type == :constant_velocity
            vx[i,j] = vx_bc_lower
        end
    end
    return nothing
end

end # module