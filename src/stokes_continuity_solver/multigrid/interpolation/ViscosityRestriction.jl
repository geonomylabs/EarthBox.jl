"""
    ViscosityRestriction

This module provides functionality for viscosity restriction operations in 3D multigrid methods.
It interpolates shear (etaxy, etaxz, etayz) and normal (etan) viscosity from finer (n) to 
coarser (n+1) level using bilinear interpolation and arithmetic averaging.

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
# Basic (density, shear viscosity - etas) nodes are shown with +
# normal viscosity (etan) is defined in pressure nodes P
# Ghost nodes shown outside the basic grid are used for boundary conditions
"""
module ViscosityRestriction

import EarthBox.Arrays.ArrayTypes.ScalarArray3D: grid_array3D
import EarthBox.Arrays.ArrayTypes.ScalarArray2D: grid_array2D
import ..TriLinearInterpolation: add_to_numerator_and_denominator!
import ..BiLinearInterpolation: add_to_numerator_and_denominator_2d!
import ..LevelManager: LevelData, LevelData2d
import ..ArrayStats

"""
    viscosity_restriction3d!(n, level_vector)

Performs restriction operation: interpolates shear (etaxy, etaxz, etayz) and normal (etan) 
viscosity from finer (n) to coarser (n+1) level using bilinear interpolation and arithmetic 
averaging.

# Arguments
- `n::Int`: Current grid level
- `level_vector::Vector{LevelData}`: Vector of level data

# Returns
- `etaxyc::Array{Float64,3}`: Interpolated shear viscosity etaxy on coarser grid
- `etaxzc::Array{Float64,3}`: Interpolated shear viscosity etaxz on coarser grid
- `etayzc::Array{Float64,3}`: Interpolated shear viscosity etayz on coarser grid
- `etanc::Array{Float64,3}`: Interpolated normal viscosity etan on coarser grid
"""
function viscosity_restriction3d!(
    n::Int, 
    level_vector::Vector{LevelData}
)::Nothing
    (
        level_vector[n+1].etaxy.array,
        level_vector[n+1].etaxz.array, 
        level_vector[n+1].etayz.array,
        level_vector[n+1].etan.array
    ) = interpolate_viscosity_to_coarser_level!(n, level_vector)
    return nothing
end

function viscosity_restriction2d!(
    n::Int, 
    level_vector::Vector{LevelData2d}
)::Nothing
    (
        level_vector[n+1].etas.array,
        level_vector[n+1].etan.array
    ) = interpolate_viscosity_to_coarser_level!(n, level_vector)
    return nothing
end

function interpolate_viscosity_to_coarser_level!(
    n::Int,
    level_vector::Vector{LevelData}
)::Tuple{Array{Float64,3}, Array{Float64,3}, Array{Float64,3}, Array{Float64,3}}
    # Coarser grid on which viscosity will be interpolated
    gridc = level_vector[n+1].grid
    xnumc = gridc.parameters.geometry.xnum.value
    ynumc = gridc.parameters.geometry.ynum.value
    znumc = gridc.parameters.geometry.znum.value
    # Creating arrays for the coarser level that act initially as numerators of 
    # trilinear interpolation equation
    # Shear and normal viscosity
    etaxyc = grid_array3D(ynumc, xnumc, znumc, Val(:shearxy))
    etaxzc = grid_array3D(ynumc, xnumc, znumc, Val(:shearxz))
    etayzc = grid_array3D(ynumc, xnumc, znumc, Val(:shearyz))
    etanc = grid_array3D(ynumc, xnumc, znumc, Val(:pressure))
    # Interpolation weights (denominators of trilinear interpolation equation)
    wtxy = copy(etaxyc)
    wtxz = copy(etaxzc)
    wtyz = copy(etayzc)
    wtn = copy(etanc)
    calculate_numerator_and_denominator_for_trilinear_interpolation!(
        n, level_vector, etaxyc, etaxzc, etayzc, etanc, wtxy, wtxz, wtyz, wtn)
    calculate_viscosity_on_coarser_level!(
        n, level_vector, etaxyc, etaxzc, etayzc, etanc, wtxy, wtxz, wtyz, wtn)
    return etaxyc, etaxzc, etayzc, etanc
end

function interpolate_viscosity_to_coarser_level!(
    n::Int,
    level_vector::Vector{LevelData2d}
)::Tuple{Array{Float64,2}, Array{Float64,2}}
    # Coarser grid on which viscosity will be interpolated
    gridc = level_vector[n+1].grid
    xnumc = gridc.parameters.geometry.xnum.value
    ynumc = gridc.parameters.geometry.ynum.value
    # Creating arrays for the coarser level that act initially as numerators of 
    # bilinear interpolation equation
    # Shear and normal viscosity
    etas = grid_array2D(ynumc, xnumc, Val(:basic))
    etan = grid_array2D(ynumc, xnumc, Val(:pressure))
    # Interpolation weights (denominators of bilinear interpolation equation)
    wts = copy(etas)
    wtn = copy(etan)
    calculate_numerator_and_denominator_for_bilinear_interpolation!(
        n, level_vector, etas, etan, wts, wtn)
    calculate_viscosity_on_coarser_level!(
        n, level_vector, etas, etan, wts, wtn)
    return etas, etan
end

function calculate_numerator_and_denominator_for_trilinear_interpolation!(
    n::Int, 
    level_vector::Vector{LevelData},
    etaxyc::Array{Float64,3},
    etaxzc::Array{Float64,3},
    etayzc::Array{Float64,3},
    etanc::Array{Float64,3},
    wtxy::Array{Float64,3},
    wtxz::Array{Float64,3},
    wtyz::Array{Float64,3},
    wtn::Array{Float64,3}
)::Nothing
    gridf = level_vector[n].grid # finer grid (current grid)

    xnumf = gridf.parameters.geometry.xnum.value
    ynumf = gridf.parameters.geometry.ynum.value
    znumf = gridf.parameters.geometry.znum.value

    etaxy = level_vector[n].etaxy.array
    etaxz = level_vector[n].etaxz.array
    etayz = level_vector[n].etayz.array
    etan = level_vector[n].etan.array

    fine_to_coarse_mapping = level_vector[n].fine_to_coarse_mapping
    etaxy_map = fine_to_coarse_mapping.etaxy_map
    etaxz_map = fine_to_coarse_mapping.etaxz_map
    etayz_map = fine_to_coarse_mapping.etayz_map
    etan_map = fine_to_coarse_mapping.pr_map

    # Calculate numerator and denominator of trilinear interpolation equation via summation
    for i in 1:ynumf
        for j in 1:xnumf
            for k in 1:znumf
                # Shear viscosity (etaxy)
                if k < znumf
                    add_to_numerator_and_denominator!(i, j, k, etaxy_map, etaxyc, wtxy, etaxy)
                end
                # Shear viscosity (etaxz)
                if i < ynumf
                    add_to_numerator_and_denominator!(i, j, k, etaxz_map, etaxzc, wtxz, etaxz)
                end
                # Shear viscosity (etayz)
                if j < xnumf
                    add_to_numerator_and_denominator!(i, j, k, etayz_map, etayzc, wtyz, etayz)
                end
                # Normal viscosity in P nodes
                if i < ynumf && j < xnumf && k < znumf
                    add_to_numerator_and_denominator!(i, j, k, etan_map, etanc, wtn, etan)
                end
            end
        end            
    end
    return nothing
end

function calculate_numerator_and_denominator_for_bilinear_interpolation!(
    n::Int, 
    level_vector::Vector{LevelData2d},
    etasc::Array{Float64,2},
    etanc::Array{Float64,2},
    wts::Array{Float64,2},
    wtn::Array{Float64,2}
)::Nothing
    gridf = level_vector[n].grid # finer grid (current grid)

    xnumf = gridf.parameters.geometry.xnum.value
    ynumf = gridf.parameters.geometry.ynum.value

    etas = level_vector[n].etas.array
    etan = level_vector[n].etan.array

    fine_to_coarse_mapping = level_vector[n].fine_to_coarse_mapping
    etas_map = fine_to_coarse_mapping.etas_map
    etan_map = fine_to_coarse_mapping.pr_map

    # Calculate numerator and denominator of bilinear interpolation equation via summation
    for i in 1:ynumf
        for j in 1:xnumf
            # Shear viscosity (etas)
            add_to_numerator_and_denominator_2d!(i, j, etas_map, etasc, wts, etas)
            # Normal viscosity (etan)
            if i < ynumf && j < xnumf
                add_to_numerator_and_denominator_2d!(i, j, etan_map, etanc, wtn, etan)
            end
        end
    end
    return nothing
end

function calculate_viscosity_on_coarser_level!(
    n::Int,
    level_vector::Vector{LevelData},
    etaxyc::Array{Float64,3},
    etaxzc::Array{Float64,3},
    etayzc::Array{Float64,3},
    etanc::Array{Float64,3},
    wtxy::Array{Float64,3},
    wtxz::Array{Float64,3},
    wtyz::Array{Float64,3},
    wtn::Array{Float64,3}
)::Nothing
    gridc = level_vector[n+1].grid # coarser grid on which viscosity will be interpolated
    xnumc = gridc.parameters.geometry.xnum.value
    ynumc = gridc.parameters.geometry.ynum.value
    znumc = gridc.parameters.geometry.znum.value
    # Recomputing viscosities for the coarser level (n+1)
    for ic in 1:ynumc
        for jc in 1:xnumc
            for kc in 1:znumc
                # Shear viscosity (etaxy)
                if kc < znumc
                    # Check for non-zero weights 
                    if wtxy[ic,jc,kc] != 0
                        etaxyc[ic,jc,kc] = etaxyc[ic,jc,kc] / wtxy[ic,jc,kc]
                    end
                end
                # Shear viscosity (etaxz)
                if ic < ynumc
                    # Check for non-zero weights 
                    if wtxz[ic,jc,kc] != 0
                        etaxzc[ic,jc,kc] = etaxzc[ic,jc,kc] / wtxz[ic,jc,kc]
                    end
                end
                # Shear viscosity (etayz)
                if jc < xnumc
                    # Check for non-zero weights 
                    if wtyz[ic,jc,kc] != 0
                        etayzc[ic,jc,kc] = etayzc[ic,jc,kc] / wtyz[ic,jc,kc]
                    end
                end
                # Normal viscosity (etan)
                if ic < ynumc && jc < xnumc && kc < znumc
                    # Check for non-zero weights 
                    if wtn[ic,jc,kc] != 0
                        etanc[ic,jc,kc] = etanc[ic,jc,kc] / wtn[ic,jc,kc]
                    end
                end
            end
        end            
    end
    return nothing
end

function calculate_viscosity_on_coarser_level!(
    n::Int,
    level_vector::Vector{LevelData2d},
    etasc::Array{Float64,2},
    etanc::Array{Float64,2},
    wts::Array{Float64,2},
    wtn::Array{Float64,2}
)::Nothing
    gridc = level_vector[n+1].grid # coarser grid on which viscosity will be interpolated
    xnumc = gridc.parameters.geometry.xnum.value
    ynumc = gridc.parameters.geometry.ynum.value
    for jc in 1:xnumc
        for ic in 1:ynumc
            if wts[ic,jc] != 0
                etasc[ic,jc] = etasc[ic,jc] / wts[ic,jc]
            end
            if ic < ynumc && jc < xnumc
                if wtn[ic,jc] != 0
                    etanc[ic,jc] = etanc[ic,jc] / wtn[ic,jc]
                end
            end
        end            
    end
    return nothing
end

end # module 