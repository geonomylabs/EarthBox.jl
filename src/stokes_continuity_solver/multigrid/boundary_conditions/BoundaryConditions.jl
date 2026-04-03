module BoundaryConditions

include("core/VxBoundaryConditions.jl")
include("core/VyBoundaryConditions.jl")
include("core/VzBoundaryConditions.jl")

import .VxBoundaryConditions: set_vx_boundary_conditions2d!
import .VyBoundaryConditions: set_vy_boundary_conditions2d!
import .VxBoundaryConditions: set_vx_boundary_conditions3d!
import .VyBoundaryConditions: set_vy_boundary_conditions3d!
import .VzBoundaryConditions: set_vz_boundary_conditions3d!
import ..LevelManager: LevelData, LevelData2d
import ..Domain: on_vx_boundary2d, on_vy_boundary2d
import ..Domain: on_vx_boundary3d, on_vy_boundary3d, on_vz_boundary3d
import ..Domain: get_ynum_vx_grid, get_znum_vx_grid
import ..Domain: get_xnum_vy_grid, get_znum_vy_grid
import ..Domain: get_ynum_vx_grid, get_znum_vx_grid
import ..Domain: get_xnum_vy_grid, get_znum_vy_grid
import ..Domain: get_xnum_vz_grid, get_ynum_vz_grid
import ..ArrayStats

function set_boundary_conditions2d!(
    level_data::LevelData2d
)::Nothing
    xnum = level_data.grid.parameters.geometry.xnum.value
    ynum = level_data.grid.parameters.geometry.ynum.value
    ynum_vx = get_ynum_vx_grid(ynum)
    xnum_vy = get_xnum_vy_grid(xnum)

    dims = (
        xnum=xnum, 
        ynum=ynum, 
        ynum_vx=ynum_vx, 
        xnum_vy=xnum_vy
        )

    Threads.@threads for j = 1:xnum+1
        for i = 1:ynum+1
            if j < xnum+1 # vx is not defined beyond xnum
                if on_vx_boundary2d(i, j, ynum, xnum)
                    set_vx_boundary_conditions2d!(i,j, dims, level_data)
                end
            end
            if i < ynum+1 # vy is not defined beyond ynum
                if on_vy_boundary2d(i, j, ynum, xnum)
                    set_vy_boundary_conditions2d!(i, j, dims, level_data)
                end
            end
        end
    end
    return nothing
end

function set_boundary_conditions2d_opt!(level_data::LevelData2d)::Nothing
    xnum = level_data.grid.parameters.geometry.xnum.value
    ynum = level_data.grid.parameters.geometry.ynum.value
    ynum_vx = get_ynum_vx_grid(ynum)
    xnum_vy = get_xnum_vy_grid(xnum)

    dims = (
        xnum=xnum, 
        ynum=ynum, 
        ynum_vx=ynum_vx, 
        xnum_vy=xnum_vy
        )

    # Set Vx boundary conditions - only on actual boundaries
    # Top boundary (i=1) and bottom boundary (i=ynum_vx) for all j
    @inbounds for j = 1:xnum
        set_vx_boundary_conditions2d!(1, j, dims, level_data)        # Top
        set_vx_boundary_conditions2d!(ynum_vx, j, dims, level_data)  # Bottom
    end
    # Left boundary (j=1) and right boundary (j=xnum) for internal i
    @inbounds for i = 2:ynum_vx-1
        set_vx_boundary_conditions2d!(i, 1, dims, level_data)     # Left
        set_vx_boundary_conditions2d!(i, xnum, dims, level_data)  # Right
    end

    # Set Vy boundary conditions - only on actual boundaries
    # Top boundary (i=1) and bottom boundary (i=ynum) for all j
    @inbounds for j = 1:xnum_vy
        set_vy_boundary_conditions2d!(1, j, dims, level_data)     # Top
        set_vy_boundary_conditions2d!(ynum, j, dims, level_data)  # Bottom
    end
    # Left boundary (j=1) and right boundary (j=xnum_vy) for internal i
    @inbounds for i = 2:ynum-1
        set_vy_boundary_conditions2d!(i, 1, dims, level_data)        # Left
        set_vy_boundary_conditions2d!(i, xnum_vy, dims, level_data)  # Right
    end
end

function set_boundary_conditions3d!(level_data::LevelData)::Nothing
    xnum = level_data.grid.parameters.geometry.xnum.value
    ynum = level_data.grid.parameters.geometry.ynum.value
    znum = level_data.grid.parameters.geometry.znum.value

    ynum_vx = get_ynum_vx_grid(ynum)
    znum_vx = get_znum_vx_grid(znum)
    xnum_vy = get_xnum_vy_grid(xnum)
    znum_vy = get_znum_vy_grid(znum)
    xnum_vz = get_xnum_vz_grid(xnum)
    ynum_vz = get_ynum_vz_grid(ynum)

    dims = (
        xnum=xnum, ynum=ynum, znum=znum,
        ynum_vx=ynum_vx, znum_vx=znum_vx,
        xnum_vy=xnum_vy, znum_vy=znum_vy,
        xnum_vz=xnum_vz, ynum_vz=ynum_vz
        )

    # Vx boundary: valid range i in 1:ynum_vx, j in 1:xnum, k in 1:znum_vx
    # Visit each boundary cell exactly once via 6 non-overlapping face loops.
    @inbounds begin
    for j = 1:xnum
        for i = 1:ynum_vx
            set_vx_boundary_conditions3d!(i, j, 1, dims, level_data)
            set_vx_boundary_conditions3d!(i, j, znum_vx, dims, level_data)
        end
    end
    for k = 2:znum_vx-1
        for j = 1:xnum
            set_vx_boundary_conditions3d!(1, j, k, dims, level_data)
            set_vx_boundary_conditions3d!(ynum_vx, j, k, dims, level_data)
        end
    end
    for k = 2:znum_vx-1
        for i = 2:ynum_vx-1
            set_vx_boundary_conditions3d!(i, 1, k, dims, level_data)
            set_vx_boundary_conditions3d!(i, xnum, k, dims, level_data)
        end
    end

    # Vy boundary: valid range i in 1:ynum, j in 1:xnum_vy, k in 1:znum_vy
    for j = 1:xnum_vy
        for i = 1:ynum
            set_vy_boundary_conditions3d!(i, j, 1, dims, level_data)
            set_vy_boundary_conditions3d!(i, j, znum_vy, dims, level_data)
        end
    end
    for k = 2:znum_vy-1
        for j = 1:xnum_vy
            set_vy_boundary_conditions3d!(1, j, k, dims, level_data)
            set_vy_boundary_conditions3d!(ynum, j, k, dims, level_data)
        end
    end
    for k = 2:znum_vy-1
        for i = 2:ynum-1
            set_vy_boundary_conditions3d!(i, 1, k, dims, level_data)
            set_vy_boundary_conditions3d!(i, xnum_vy, k, dims, level_data)
        end
    end

    # Vz boundary: valid range i in 1:ynum_vz, j in 1:xnum_vz, k in 1:znum
    for j = 1:xnum_vz
        for i = 1:ynum_vz
            set_vz_boundary_conditions3d!(i, j, 1, dims, level_data)
            set_vz_boundary_conditions3d!(i, j, znum, dims, level_data)
        end
    end
    for k = 2:znum-1
        for j = 1:xnum_vz
            set_vz_boundary_conditions3d!(1, j, k, dims, level_data)
            set_vz_boundary_conditions3d!(ynum_vz, j, k, dims, level_data)
        end
    end
    for k = 2:znum-1
        for i = 2:ynum_vz-1
            set_vz_boundary_conditions3d!(i, 1, k, dims, level_data)
            set_vz_boundary_conditions3d!(i, xnum_vz, k, dims, level_data)
        end
    end
    end # @inbounds
    return nothing
end

end