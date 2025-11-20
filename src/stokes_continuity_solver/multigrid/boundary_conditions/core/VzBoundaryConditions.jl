module VzBoundaryConditions

import ...LevelManager: LevelData
import ...Domain: get_xnum_vz_grid, get_ynum_vz_grid

function set_vz_boundary_conditions3d!(
    i::Int64,
    j::Int64,
    k::Int64,
    dims::NamedTuple,
    level_data::LevelData;
    bc_type::Symbol = :free_slip,
    vz_bc_left::Float64 = 0.0,
    vz_bc_right::Float64 = 0.0,
    vz_bc_upper::Float64 = 0.0,
    vz_bc_lower::Float64 = 0.0,
    vz_bc_back::Float64 = 0.0,
    vz_bc_front::Float64 = 0.0
)::Nothing
    vz = level_data.vz.array
    znum = dims.znum
    ynum_vz = dims.ynum_vz
    xnum_vz = dims.xnum_vz
    ynum_vz = dims.ynum_vz
    
    vz[i,j,k] = 0.0
    # Back Boundary
    if k == 1
        if bc_type == :constant_velocity
            vz[i,j,k] = vz_bc_back
        end
    end
    # Front Boundary
    if k == znum
        if bc_type == :constant_velocity
            vz[i,j,k] = vz_bc_front
        end
    end
    # Upper Boundary
    if i == 1
        # Free slip
        if bc_type == :free_slip
            vz[i,j,k] = vz[i+1,j,k]
        # No slip
        elseif bc_type == :no_slip
            vz[i,j,k] = -vz[i+1,j,k]
        elseif bc_type == :constant_velocity
            vz[i,j,k] = vz_bc_upper
        end
    end
    # Lower boundary
    if i == ynum_vz
        # Free slip
        if bc_type == :free_slip
            vz[i,j,k] = vz[i-1,j,k]
        # No slip
        elseif bc_type == :no_slip
            vz[i,j,k] = -vz[i-1,j,k]
        elseif bc_type == :constant_velocity
            vz[i,j,k] = vz_bc_lower
        end
    end
    # Left boundary
    if j == 1
        # Free slip
        if bc_type == :free_slip
            vz[i,j,k] = vz[i,j+1,k]
        # No slip
        elseif bc_type == :no_slip
            vz[i,j,k] = -vz[i,j+1,k]
        elseif bc_type == :constant_velocity
            vz[i,j,k] = vz_bc_left
        end
    end
    # Right boundary
    if j == xnum_vz
        # Free slip
        if bc_type == :free_slip
            vz[i,j,k] = vz[i,j-1,k]
        # No slip
        elseif bc_type == :no_slip
            vz[i,j,k] = -vz[i,j-1,k]
        elseif bc_type == :constant_velocity
            vz[i,j,k] = vz_bc_right
        end
    end
    return nothing
end

end # module