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

    Threads.@threads for k = 1:znum+1
        for j = 1:xnum+1
            for i = 1:ynum+1
                if j < xnum+1 # vx is not defined beyond xnum
                    if on_vx_boundary3d(i, j, k, ynum, xnum, znum)
                        set_vx_boundary_conditions3d!(i,j,k, dims, level_data)
                    end
                end
                if i < ynum+1 # vy is not defined beyond ynum
                    if on_vy_boundary3d(i, j, k, ynum, xnum, znum)
                        set_vy_boundary_conditions3d!(i, j, k, dims, level_data)
                    end
                end
                if k < znum+1 # vz is not defined beyond znum
                    if on_vz_boundary3d(i, j, k, ynum, xnum, znum)
                        set_vz_boundary_conditions3d!(i, j, k, dims, level_data)
                    end
                end
            end
        end
    end
    return nothing
end


end