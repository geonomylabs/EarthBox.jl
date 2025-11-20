"""
    Restriction

This module provides the function `stokes_continuity3d_viscous_restriction()` which performs 
restriction operations for 3D Stokes continuity equations with viscous terms.

The function interpolates residuals (resx, resy, resz, resc) from finer (n) to coarser (n+1) 
level using trilinear interpolation and produces right parts (RX, RY, RZ, RC) for the coarser level.
Normal viscosity on finer (etanf) and coarser (etanc) levels is used for rescaling continuity residuals.

Resolution (xnum, ynum, znum) and steps (xstp, ystp, zstp) for both levels are used for 
organizing the interpolation.

# Staggered Grid for Multigrid
# 
#     vx       vx       vx    
#
# vy  +---vy---+---vy---+   vy
#     |        |        |
#     vx   P   vx   P   vx    
#     |        |        |
# vy  +---vy---+---vy---+   vy
#     |        |        |
#     vx   P   vx   P   vx    
#     |        |        |
# vy  +---vy---+---vy---+   vy
#
#     vx       vx       vx    
# 
# Lines show basic grid
# Basic (density) nodes are shown with +
# Ghost nodes shown outside the basic grid
# are used for boundary conditions
"""
module Restriction

import EarthBox.Arrays.ArrayTypes.ScalarArray3D: grid_array3D
import EarthBox.Arrays.ArrayTypes.ScalarArray2D: grid_array2D
import ..TriLinearInterpolation: add_to_numerator_and_denominator!
import ..BiLinearInterpolation: add_to_numerator_and_denominator_2d!
import ..LevelManager: LevelData, LevelData2d

"""
    restrict_stokes3d_residuals!(n, level_vector, resx, resy, resz, resc)

Perform restriction operation for 3D Stokes continuity equations with viscous terms.

# Arguments
- `n::Int`: Current grid level
- `level_vector::Vector{LevelData}`: Vector of level data
- `resx::Array{Float64,3}`: x-Stokes equation residuals on finer level
- `resy::Array{Float64,3}`: y-Stokes equation residuals on finer level
- `resz::Array{Float64,3}`: z-Stokes equation residuals on finer level
- `resc::Array{Float64,3}`: Continuity equation residuals on finer level

"""
function restrict_stokes3d_residuals!(
    n::Int64,
    level_vector::Vector{LevelData},
    ΔRx_fine::Array{Float64,3}, 
    ΔRy_fine::Array{Float64,3}, 
    ΔRz_fine::Array{Float64,3}, 
    ΔRc_fine::Array{Float64,3},
)::Nothing
    (
        level_vector[n+1].RX.array,
        level_vector[n+1].RY.array,
        level_vector[n+1].RZ.array,
        level_vector[n+1].RC.array
    ) = interpolate_residuals_to_coarser_level!(
        n, level_vector, ΔRx_fine, ΔRy_fine, ΔRz_fine, ΔRc_fine)
    return nothing
end

function interpolate_residuals_to_coarser_level!(
    n::Int64,
    level_vector::Vector{LevelData},
    ΔRx_fine::Array{Float64,3}, 
    ΔRy_fine::Array{Float64,3}, 
    ΔRz_fine::Array{Float64,3}, 
    ΔRc_fine::Array{Float64,3},
)::Tuple{Array{Float64,3}, Array{Float64,3}, Array{Float64,3}, Array{Float64,3}}
    gridc = level_vector[n+1].grid
    xnumc = gridc.parameters.geometry.xnum.value
    ynumc = gridc.parameters.geometry.ynum.value
    znumc = gridc.parameters.geometry.znum.value
    # Creating arrays for the coarser level
    # Right parts
    ΔRx_coarse = grid_array3D(ynumc, xnumc, znumc, Val(:vx))
    ΔRy_coarse = grid_array3D(ynumc, xnumc, znumc, Val(:vy))
    ΔRz_coarse = grid_array3D(ynumc, xnumc, znumc, Val(:vz))
    ΔRc_coarse = grid_array3D(ynumc, xnumc, znumc, Val(:pressure))
    # Interpolation weights
    wtx = grid_array3D(ynumc, xnumc, znumc, Val(:vx))
    wty = grid_array3D(ynumc, xnumc, znumc, Val(:vy))
    wtz = grid_array3D(ynumc, xnumc, znumc, Val(:vz))
    wtc = grid_array3D(ynumc, xnumc, znumc, Val(:pressure))
    calculate_numerator_and_denominator_for_trilinear_interpolation!(
        n, level_vector,
        ΔRx_fine, ΔRy_fine, ΔRz_fine, ΔRc_fine, 
        ΔRx_coarse, ΔRy_coarse, ΔRz_coarse, ΔRc_coarse, 
        wtx, wty, wtz, wtc
        )
    calculate_residuals_on_coarser_level!(
        n, level_vector, 
        ΔRx_coarse, ΔRy_coarse, ΔRz_coarse, ΔRc_coarse, 
        wtx, wty, wtz, wtc
        )
    return ΔRx_coarse, ΔRy_coarse, ΔRz_coarse, ΔRc_coarse
end

function calculate_numerator_and_denominator_for_trilinear_interpolation!(
    n::Int64,
    level_vector::Vector{LevelData},
    ΔRx_fine::Array{Float64,3}, 
    ΔRy_fine::Array{Float64,3}, 
    ΔRz_fine::Array{Float64,3}, 
    ΔRc_fine::Array{Float64,3},
    ΔRx_coarse::Array{Float64,3},
    ΔRy_coarse::Array{Float64,3},
    ΔRz_coarse::Array{Float64,3},
    ΔRc_coarse::Array{Float64,3},
    wtx::Array{Float64,3},
    wty::Array{Float64,3},
    wtz::Array{Float64,3},
    wtc::Array{Float64,3},
)::Nothing
    gridf = level_vector[n].grid

    xnumf = gridf.parameters.geometry.xnum.value
    ynumf = gridf.parameters.geometry.ynum.value
    znumf = gridf.parameters.geometry.znum.value

    etanf = level_vector[n].etan.array
    etan_resc = etanf .* ΔRc_fine

    fine_to_coarse_mapping = level_vector[n].fine_to_coarse_mapping
    vx_map = fine_to_coarse_mapping.vx_map
    vy_map = fine_to_coarse_mapping.vy_map
    vz_map = fine_to_coarse_mapping.vz_map
    pr_map = fine_to_coarse_mapping.pr_map
    
    # Interpolating residuals from finer level nodes
    # Cycle of node of finer (n) level
    for k in 2:znumf
        for j in 2:xnumf
            for i in 2:ynumf
                # x-Stokes equation residual
                if j < xnumf
                    add_to_numerator_and_denominator!(i, j, k, vx_map, ΔRx_coarse, wtx, ΔRx_fine)
                end
                # y-Stokes equation residual
                if i < ynumf
                    add_to_numerator_and_denominator!(i, j, k, vy_map, ΔRy_coarse, wty, ΔRy_fine)
                end
                # z-Stokes equation residual
                if k < znumf
                    add_to_numerator_and_denominator!(i, j, k, vz_map, ΔRz_coarse, wtz, ΔRz_fine)
                end
                # Continuity equation residual
                if i < ynumf && j < xnumf && k < znumf
                    add_to_numerator_and_denominator!(i, j, k, pr_map, ΔRc_coarse, wtc, etan_resc)
                end
            
            end
        end            
    end
    return nothing
end

function calculate_residuals_on_coarser_level!(
    n::Int64,
    level_vector::Vector{LevelData},
    ΔRx_coarse::Array{Float64,3},
    ΔRy_coarse::Array{Float64,3},
    ΔRz_coarse::Array{Float64,3},
    ΔRc_coarse::Array{Float64,3},
    wtx::Array{Float64,3},
    wty::Array{Float64,3},
    wtz::Array{Float64,3},
    wtc::Array{Float64,3},
)::Nothing
    gridc = level_vector[n+1].grid
    xnumc = gridc.parameters.geometry.xnum.value
    ynumc = gridc.parameters.geometry.ynum.value
    znumc = gridc.parameters.geometry.znum.value
    etanc = level_vector[n+1].etan.array
    # Recomputing right parts (RX, RY, RZ, RC)
    # for the coarser level (n+1)
    for kc in 1:znumc+1
        for jc in 1:xnumc+1
            for ic in 1:ynumc+1
                # x-Stokes
                if jc < xnumc+1
                    # Check for non-zero weights 
                    if wtx[ic,jc,kc] != 0 && ic > 1 && ic < ynumc+1 && 
                       jc > 1 && jc < xnumc && kc > 1 && kc < znumc+1
                        ΔRx_coarse[ic,jc,kc] = ΔRx_coarse[ic,jc,kc]/wtx[ic,jc,kc]
                    else
                        ΔRx_coarse[ic,jc,kc] = 0
                    end
                end
                # y-Stokes
                if ic < ynumc+1
                    # Check for non-zero weights 
                    if wty[ic,jc,kc] != 0 && ic > 1 && ic < ynumc && 
                       jc > 1 && jc < xnumc+1 && kc > 1 && kc < znumc+1
                        ΔRy_coarse[ic,jc,kc] = ΔRy_coarse[ic,jc,kc]/wty[ic,jc,kc]
                    else
                        ΔRy_coarse[ic,jc,kc] = 0
                    end
                end
                # z-Stokes
                if kc < znumc+1
                    # Check for non-zero weights 
                    if wtz[ic,jc,kc] != 0 && ic > 1 && ic < ynumc+1 && 
                       jc > 1 && jc < xnumc+1 && kc > 1 && kc < znumc
                        ΔRz_coarse[ic,jc,kc] = ΔRz_coarse[ic,jc,kc]/wtz[ic,jc,kc]
                    else
                        ΔRz_coarse[ic,jc,kc] = 0
                    end
                end
                # Continuity
                if ic < ynumc && jc < xnumc && kc < znumc
                    if wtc[ic,jc,kc] != 0
                        ΔRc_coarse[ic,jc,kc] = ΔRc_coarse[ic,jc,kc]/wtc[ic,jc,kc]/etanc[ic,jc,kc]
                    else
                        ΔRc_coarse[ic,jc,kc] = 0
                    end
                end
            
            end
        end            
    end
    return nothing
end

function restrict_stokes2d_residuals!(
    n::Int64,
    level_vector::Vector{LevelData2d},
    ΔRx_fine::Array{Float64,2}, 
    ΔRy_fine::Array{Float64,2}, 
    ΔRc_fine::Array{Float64,2},
)::Nothing
    (
        level_vector[n+1].RX.array,
        level_vector[n+1].RY.array,
        level_vector[n+1].RC.array
    ) = interpolate_residuals_to_coarser_level_2d!(n, level_vector, ΔRx_fine, ΔRy_fine, ΔRc_fine)
    return nothing
end

function interpolate_residuals_to_coarser_level_2d!(
    n::Int64,
    level_vector::Vector{LevelData2d},
    ΔRx_fine::Array{Float64,2}, 
    ΔRy_fine::Array{Float64,2}, 
    ΔRc_fine::Array{Float64,2},
)::Tuple{Array{Float64,2}, Array{Float64,2}, Array{Float64,2}}
    gridc = level_vector[n+1].grid
    xnumc = gridc.parameters.geometry.xnum.value
    ynumc = gridc.parameters.geometry.ynum.value
    # Creating arrays for the coarser level
    # Right parts
    ΔRx_coarse = grid_array2D(ynumc, xnumc, Val(:vx))
    ΔRy_coarse = grid_array2D(ynumc, xnumc, Val(:vy))
    ΔRc_coarse = grid_array2D(ynumc, xnumc, Val(:pressure))
    # Interpolation weights
    wtx = grid_array2D(ynumc, xnumc, Val(:vx))
    wty = grid_array2D(ynumc, xnumc, Val(:vy))
    wtc = grid_array2D(ynumc, xnumc, Val(:pressure))
    calculate_numerator_and_denominator_for_bilinear_interpolation!(
        n, level_vector,
        ΔRx_fine, ΔRy_fine, ΔRc_fine, 
        ΔRx_coarse, ΔRy_coarse, ΔRc_coarse, 
        wtx, wty, wtc
        )
    calculate_residuals_on_coarser_level_2d!(
        n, level_vector, 
        ΔRx_coarse, ΔRy_coarse, ΔRc_coarse, 
        wtx, wty, wtc
        )
    return ΔRx_coarse, ΔRy_coarse, ΔRc_coarse
end

function calculate_numerator_and_denominator_for_bilinear_interpolation!(
    n::Int64,
    level_vector::Vector{LevelData2d},
    ΔRx_fine::Array{Float64,2}, 
    ΔRy_fine::Array{Float64,2}, 
    ΔRc_fine::Array{Float64,2},
    ΔRx_coarse::Array{Float64,2},
    ΔRy_coarse::Array{Float64,2},
    ΔRc_coarse::Array{Float64,2},
    wtx::Array{Float64,2},
    wty::Array{Float64,2},
    wtc::Array{Float64,2},
)::Nothing
    gridf = level_vector[n].grid

    xnumf = gridf.parameters.geometry.xnum.value
    ynumf = gridf.parameters.geometry.ynum.value

    etanf = level_vector[n].etan.array
    etan_resc = etanf .* ΔRc_fine

    fine_to_coarse_mapping = level_vector[n].fine_to_coarse_mapping
    vx_map = fine_to_coarse_mapping.vx_map
    vy_map = fine_to_coarse_mapping.vy_map
    pr_map = fine_to_coarse_mapping.pr_map
    
    for j in 2:xnumf
        for i in 2:ynumf
            if j < xnumf
                add_to_numerator_and_denominator_2d!(i, j, vx_map, ΔRx_coarse, wtx, ΔRx_fine)
            end
            if i < ynumf
                add_to_numerator_and_denominator_2d!(i, j, vy_map, ΔRy_coarse, wty, ΔRy_fine)
            end
            if i < ynumf && j < xnumf
                add_to_numerator_and_denominator_2d!(i, j, pr_map, ΔRc_coarse, wtc, etan_resc)
            end
        end            
    end
    return nothing
end

function calculate_residuals_on_coarser_level_2d!(
    n::Int64,
    level_vector::Vector{LevelData2d},
    ΔRx_coarse::Array{Float64,2},
    ΔRy_coarse::Array{Float64,2},
    ΔRc_coarse::Array{Float64,2},
    wtx::Array{Float64,2},
    wty::Array{Float64,2},
    wtc::Array{Float64,2},
)::Nothing
    gridc = level_vector[n+1].grid
    xnumc = gridc.parameters.geometry.xnum.value
    ynumc = gridc.parameters.geometry.ynum.value
    etanc = level_vector[n+1].etan.array
    for jc in 1:xnumc+1
        for ic in 1:ynumc+1
            if jc < xnumc+1
                if wtx[ic,jc] != 0 && ic > 1 && ic < ynumc+1 && jc > 1 && jc < xnumc
                    ΔRx_coarse[ic,jc] = ΔRx_coarse[ic,jc]/wtx[ic,jc]
                else
                    ΔRx_coarse[ic,jc] = 0
                end
            end
            if ic < ynumc+1
                if wty[ic,jc] != 0 && ic > 1 && ic < ynumc && jc > 1 && jc < xnumc+1
                    ΔRy_coarse[ic,jc] = ΔRy_coarse[ic,jc]/wty[ic,jc]
                else
                    ΔRy_coarse[ic,jc] = 0
                end
            end
            if ic < ynumc && jc < xnumc
                if wtc[ic,jc] != 0
                    ΔRc_coarse[ic,jc] = ΔRc_coarse[ic,jc]/wtc[ic,jc]/etanc[ic,jc]
                else
                    ΔRc_coarse[ic,jc] = 0
                end
            end
        end            
    end
    return nothing
end

end # module 