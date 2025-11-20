"""
    Prolongation

This module provides the prolongation operation for Stokes-Continuity 3D multigrid:
Interpolates corrections for solution (vx,vy,vz,pr) from coarser (n) to finer (n-1) level
using trilinear interpolation and produces updated solutions (dvx,dvy,dvz,dpr) for this finer level.

Resolution (xnum,ynum,znum) and steps (xstp,ystp,zstp) for both levels
is used for organizing the interpolation.

Staggered Grid for Multigrid:

     vx       vx       vx    

vy  +---vy---+---vy---+   vy
     |        |        |
     vx   P   vx   P   vx    
     |        |        |
vy  +---vy---+---vy---+   vy
     |        |        |
     vx   P   vx   P   vx    
     |        |        |
vy  +---vy---+---vy---+   vy

     vx       vx       vx    

Lines show basic grid
Basic (density) nodes are shown with +
Ghost nodes shown outside the basic grid
are used for boundary conditions
"""
module Prolongation

import EarthBox.Arrays.ArrayTypes.ScalarArray3D: grid_array3D
import EarthBox.Arrays.ArrayTypes.ScalarArray2D: grid_array2D
import ..TriLinearInterpolation: coarse2fine_trilinear_interpolation!
import ..BiLinearInterpolation: coarse2fine_bilinear_interpolation!
import ..LevelManager: LevelData, LevelData2d

"""
    prolongate_stokes3d_solution(n, level_vector)

Makes prolongation operation: Interpolates corrections for solution (vx,vy,vz,pr) 
from coarser (n) to finer (n-1) level using trilinear interpolation
and produces updated solutions (dvx,dvy,dvz,dpr) for this finer level.

# Arguments
- `n::Int`: Current grid level
- `level_vector::Vector{LevelData}`: Vector of level data

# Returns
- `dvx::Array{Float64,3}`: x-velocity corrections for finer grid
- `dvy::Array{Float64,3}`: y-velocity corrections for finer grid
- `dvz::Array{Float64,3}`: z-velocity corrections for finer grid
- `dpr::Array{Float64,3}`: pressure corrections for finer grid
"""
function prolongate_stokes3d_solution(
    n::Int, 
    level_vector::Vector{LevelData},
)::Tuple{Array{Float64,3}, Array{Float64,3}, Array{Float64,3}, Array{Float64,3}}
    gridf = level_vector[n-1].grid
    xnumf = gridf.parameters.geometry.xnum.value
    ynumf = gridf.parameters.geometry.ynum.value
    znumf = gridf.parameters.geometry.znum.value
    
    vx = level_vector[n].vx.array
    vy = level_vector[n].vy.array
    vz = level_vector[n].vz.array
    pr = level_vector[n].pr.array

    dvx = grid_array3D(ynumf, xnumf, znumf, Val(:vx))
    dvy = grid_array3D(ynumf, xnumf, znumf, Val(:vy))
    dvz = grid_array3D(ynumf, xnumf, znumf, Val(:vz))
    dpr = grid_array3D(ynumf, xnumf, znumf, Val(:pressure))

    fine_to_coarse_mapping = level_vector[n-1].fine_to_coarse_mapping
    vx_map = fine_to_coarse_mapping.vx_map
    vy_map = fine_to_coarse_mapping.vy_map
    vz_map = fine_to_coarse_mapping.vz_map
    pr_map = fine_to_coarse_mapping.pr_map
    
    for k = 1:znumf+1
        for j = 1:xnumf+1
            for i = 1:ynumf+1
                # x-Stokes correction
                if j < xnumf+1
                    coarse2fine_trilinear_interpolation!(i, j, k, vx_map, vx, dvx)
                end
                # y-Stokes correction
                if i < ynumf+1
                    coarse2fine_trilinear_interpolation!(i, j, k, vy_map, vy, dvy)
                end
                # z-Stokes correction
                if k < znumf+1
                    coarse2fine_trilinear_interpolation!(i, j, k, vz_map, vz, dvz)
                end
                # Continuity equation correction
                if i < ynumf && j < xnumf && k < znumf
                    coarse2fine_trilinear_interpolation!(i, j, k, pr_map, pr, dpr)
                end
            
            end
        end            
    end
    return dvx, dvy, dvz, dpr
end

function prolongate_stokes2d_solution(
    n::Int, 
    level_vector::Vector{LevelData2d},
)::Tuple{Array{Float64,2}, Array{Float64,2}, Array{Float64,2}}
    gridf = level_vector[n-1].grid
    xnumf = gridf.parameters.geometry.xnum.value
    ynumf = gridf.parameters.geometry.ynum.value
    
    vx = level_vector[n].vx.array
    vy = level_vector[n].vy.array
    pr = level_vector[n].pr.array

    dvx = grid_array2D(ynumf, xnumf, Val(:vx))
    dvy = grid_array2D(ynumf, xnumf, Val(:vy))
    dpr = grid_array2D(ynumf, xnumf, Val(:pressure))

    fine_to_coarse_mapping = level_vector[n-1].fine_to_coarse_mapping
    vx_map = fine_to_coarse_mapping.vx_map
    vy_map = fine_to_coarse_mapping.vy_map
    pr_map = fine_to_coarse_mapping.pr_map
    
    for j = 1:xnumf+1
        for i = 1:ynumf+1
            # x-Stokes correction
            if j < xnumf+1
                coarse2fine_bilinear_interpolation!(i, j, vx_map, vx, dvx)
            end
            # y-Stokes correction
            if i < ynumf+1
                coarse2fine_bilinear_interpolation!(i, j, vy_map, vy, dvy)
            end
            # Continuity equation correction
            if i < ynumf && j < xnumf
                coarse2fine_bilinear_interpolation!(i, j, pr_map, pr, dpr)
            end
        end            
    end
    return dvx, dvy, dpr
end

end # module 