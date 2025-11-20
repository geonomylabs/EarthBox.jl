module VyBoundaryConditions

import ...LevelManager: LevelData, LevelData2d
import ...Domain: get_xnum_vy_grid, get_znum_vy_grid

function set_vy_boundary_conditions3d!(
    i::Int64,
    j::Int64,
    k::Int64,
    dims::NamedTuple,
    level_data::LevelData;
    bc_type::Symbol = :free_slip,
    vy_bc_left::Float64 = 0.0,
    vy_bc_right::Float64 = 0.0,
    vy_bc_upper::Float64 = 0.0,
    vy_bc_lower::Float64 = 0.0,
    vy_bc_back::Float64 = 0.0,
    vy_bc_front::Float64 = 0.0
)::Nothing
    vy = level_data.vy.array
    ynum = dims.ynum
    xnum_vy = dims.xnum_vy
    znum_vy = dims.znum_vy

    vy[i,j,k] = 0.0
    # Upper Boundary
    if i == 1
        if bc_type == :constant_velocity
            vy[i,j,k] = vy_bc_upper
        end
    end
    # Lower Boundary
    if i == ynum
        if bc_type == :constant_velocity
            vy[i,j,k] = vy_bc_lower
        end
    end
    # Left boundary
    if j == 1
        # Free slip
        if bc_type == :free_slip
            vy[i,j,k] = vy[i,j+1,k]
        # No slip
        elseif bc_type == :no_slip
            vy[i,j,k] = -vy[i,j+1,k]
        elseif bc_type == :constant_velocity
            vy[i,j,k] = vy_bc_left
        end
    end
    # Right boundary
    if j == xnum_vy
        # Free slip
        if bc_type == :free_slip
            vy[i,j,k] = vy[i,j-1,k]
        # No slip
        elseif bc_type == :no_slip
            vy[i,j,k] = -vy[i,j-1,k]
        elseif bc_type == :constant_velocity
            vy[i,j,k] = vy_bc_right
        end
    end
    # Back Boundary
    if k == 1
        # Free slip
        if bc_type == :free_slip
            vy[i,j,k] = vy[i,j,k+1]
        # No slip
        elseif bc_type == :no_slip
            vy[i,j,k] = -vy[i,j,k+1]
        elseif bc_type == :constant_velocity
            vy[i,j,k] = vy_bc_back
        end
    end
    # Front boundary
    if k == znum_vy
        # Free slip
        if bc_type == :free_slip
            vy[i,j,k] = vy[i,j,k-1]
        # No slip
        elseif bc_type == :no_slip
            vy[i,j,k] = -vy[i,j,k-1]
        elseif bc_type == :constant_velocity
            vy[i,j,k] = vy_bc_front
        end
    end
    return nothing
end

function set_vy_boundary_conditions2d!(
    i::Int64,
    j::Int64,
    dims::NamedTuple,
    level_data::LevelData2d;
    bc_type::Symbol = :free_slip,
    vy_bc_left::Float64 = 0.0,
    vy_bc_right::Float64 = 0.0,
    vy_bc_upper::Float64 = 0.0,
    vy_bc_lower::Float64 = 0.0
)::Nothing
    vy = level_data.vy.array
    ynum = dims.ynum
    xnum_vy = dims.xnum_vy

    vy[i,j] = 0.0
    # Upper Boundary
    if i == 1
        if bc_type == :constant_velocity
            vy[i,j] = vy_bc_upper
        end
    end
    # Lower Boundary
    if i == ynum
        if bc_type == :constant_velocity
            vy[i,j] = vy_bc_lower
        end
    end
    # Left boundary
    if j == 1
        # Free slip
        if bc_type == :free_slip
            vy[i,j] = vy[i,j+1]
        # No slip
        elseif bc_type == :no_slip
            vy[i,j] = -vy[i,j+1]
        elseif bc_type == :constant_velocity
            vy[i,j] = vy_bc_left
        end
    end
    # Right boundary
    if j == xnum_vy
        # Free slip
        if bc_type == :free_slip
            vy[i,j] = vy[i,j-1]
        # No slip
        elseif bc_type == :no_slip
            vy[i,j] = -vy[i,j-1]
        elseif bc_type == :constant_velocity
            vy[i,j] = vy_bc_right
        end
    end
    return nothing
end

end # module