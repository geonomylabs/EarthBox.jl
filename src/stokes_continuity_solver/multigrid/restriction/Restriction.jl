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

# Match SolveStokes3dOpt / Residuals: avoid thread overhead on tiny coarse grids.
const RESTRICT_DIVIDE_K_PARALLEL_MIN_ZNUM = 16

# Parallel restriction accumulation: partition fine k-planes across threads; each thread sums into
# coarse-sized buffers (LevelData.restrict_thread_accum), then reduce into RX/weights.
const RESTRICT_NUMER_K_PARALLEL_MIN_ZNUM = 16

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
    interpolate_residuals_to_coarser_level!(
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
)::Nothing
    coarse_ld = level_vector[n+1]
    ΔRx_coarse = coarse_ld.RX.array
    ΔRy_coarse = coarse_ld.RY.array
    ΔRz_coarse = coarse_ld.RZ.array
    ΔRc_coarse = coarse_ld.RC.array
    fill!(ΔRx_coarse, 0.0)
    fill!(ΔRy_coarse, 0.0)
    fill!(ΔRz_coarse, 0.0)
    fill!(ΔRc_coarse, 0.0)
    wtx = coarse_ld.restrict_wtx
    wty = coarse_ld.restrict_wty
    wtz = coarse_ld.restrict_wtz
    wtc = coarse_ld.restrict_wtc
    fill!(wtx, 0.0)
    fill!(wty, 0.0)
    fill!(wtz, 0.0)
    fill!(wtc, 0.0)
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
    return nothing
end

function _restrict_numer_fine_k_q_ranges(n_k::Int, nt::Int)::Vector{Tuple{Int,Int}}
    base = div(n_k, nt)
    extra = mod(n_k, nt)
    offset = 0
    ranges = Vector{Tuple{Int,Int}}(undef, nt)
    @inbounds for t in 1:nt
        len = base + (t <= extra ? 1 : 0)
        q_start = offset + 1
        q_end = offset + len
        ranges[t] = (q_start, q_end)
        offset += len
    end
    return ranges
end

function _calculate_numerator_one_k_range!(
    k_lo::Int,
    k_hi::Int,
    xnumf::Int,
    ynumf::Int,
    znumf::Int,
    vx_map,
    vy_map,
    vz_map,
    pr_map,
    ΔRx_fine::Array{Float64,3},
    ΔRy_fine::Array{Float64,3},
    ΔRz_fine::Array{Float64,3},
    etan_resc::Array{Float64,3},
    ΔRx_coarse::Array{Float64,3},
    ΔRy_coarse::Array{Float64,3},
    ΔRz_coarse::Array{Float64,3},
    ΔRc_coarse::Array{Float64,3},
    wtx::Array{Float64,3},
    wty::Array{Float64,3},
    wtz::Array{Float64,3},
    wtc::Array{Float64,3},
)::Nothing
    @inbounds for k in k_lo:k_hi
        for j in 2:xnumf
            for i in 2:ynumf
                if j < xnumf
                    add_to_numerator_and_denominator!(i, j, k, vx_map, ΔRx_coarse, wtx, ΔRx_fine)
                end
                if i < ynumf
                    add_to_numerator_and_denominator!(i, j, k, vy_map, ΔRy_coarse, wty, ΔRy_fine)
                end
                if k < znumf
                    add_to_numerator_and_denominator!(i, j, k, vz_map, ΔRz_coarse, wtz, ΔRz_fine)
                end
                if i < ynumf && j < xnumf && k < znumf
                    add_to_numerator_and_denominator!(i, j, k, pr_map, ΔRc_coarse, wtc, etan_resc)
                end
            end
        end
    end
    return nothing
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
    etan_resc = level_vector[n].etan_resc_buf
    @. etan_resc = etanf * ΔRc_fine

    fine_to_coarse_mapping = level_vector[n].fine_to_coarse_mapping
    vx_map = fine_to_coarse_mapping.vx_map
    vy_map = fine_to_coarse_mapping.vy_map
    vz_map = fine_to_coarse_mapping.vz_map
    pr_map = fine_to_coarse_mapping.pr_map

    coarse_ld = level_vector[n+1]
    bufs = coarse_ld.restrict_thread_accum
    n_k = znumf - 1
    nt_req = Threads.nthreads()
    parallel_k =
        nt_req > 1 &&
        znumf >= RESTRICT_NUMER_K_PARALLEL_MIN_ZNUM &&
        bufs !== nothing &&
        n_k >= 1 &&
        length(bufs) >= nt_req

    if !parallel_k
        _calculate_numerator_one_k_range!(
            2, znumf, xnumf, ynumf, znumf,
            vx_map, vy_map, vz_map, pr_map,
            ΔRx_fine, ΔRy_fine, ΔRz_fine, etan_resc,
            ΔRx_coarse, ΔRy_coarse, ΔRz_coarse, ΔRc_coarse,
            wtx, wty, wtz, wtc,
        )
        return nothing
    end

    nt = nt_req
    ranges = _restrict_numer_fine_k_q_ranges(n_k, nt)
    Threads.@threads for t in 1:nt
        q_start, q_end = ranges[t]
        b = bufs[t]
        ΔRx_th, wtx_th, ΔRy_th, wty_th, ΔRz_th, wtz_th, ΔRc_th, wtc_th = b
        fill!(ΔRx_th, 0.0)
        fill!(wtx_th, 0.0)
        fill!(ΔRy_th, 0.0)
        fill!(wty_th, 0.0)
        fill!(ΔRz_th, 0.0)
        fill!(wtz_th, 0.0)
        fill!(ΔRc_th, 0.0)
        fill!(wtc_th, 0.0)
        if q_end < q_start
            continue
        end
        k_lo = q_start + 1
        k_hi = q_end + 1
        _calculate_numerator_one_k_range!(
            k_lo, k_hi, xnumf, ynumf, znumf,
            vx_map, vy_map, vz_map, pr_map,
            ΔRx_fine, ΔRy_fine, ΔRz_fine, etan_resc,
            ΔRx_th, ΔRy_th, ΔRz_th, ΔRc_th,
            wtx_th, wty_th, wtz_th, wtc_th,
        )
    end

    b1 = bufs[1]
    copyto!(ΔRx_coarse, b1[1])
    copyto!(wtx, b1[2])
    copyto!(ΔRy_coarse, b1[3])
    copyto!(wty, b1[4])
    copyto!(ΔRz_coarse, b1[5])
    copyto!(wtz, b1[6])
    copyto!(ΔRc_coarse, b1[7])
    copyto!(wtc, b1[8])
    @inbounds for t in 2:nt
        bt = bufs[t]
        @. ΔRx_coarse += bt[1]
        @. wtx += bt[2]
        @. ΔRy_coarse += bt[3]
        @. wty += bt[4]
        @. ΔRz_coarse += bt[5]
        @. wtz += bt[6]
        @. ΔRc_coarse += bt[7]
        @. wtc += bt[8]
    end
    return nothing
end

function _divide_restricted_residuals_on_coarser_kplane!(
    kc::Int,
    xnumc::Int,
    ynumc::Int,
    znumc::Int,
    ΔRx_coarse::Array{Float64,3},
    ΔRy_coarse::Array{Float64,3},
    ΔRz_coarse::Array{Float64,3},
    ΔRc_coarse::Array{Float64,3},
    wtx::Array{Float64,3},
    wty::Array{Float64,3},
    wtz::Array{Float64,3},
    wtc::Array{Float64,3},
    etanc::Array{Float64,3},
)::Nothing
    @inbounds for jc in 1:xnumc+1
        for ic in 1:ynumc+1
            # x-Stokes
            if jc < xnumc+1
                if wtx[ic,jc,kc] != 0 && ic > 1 && ic < ynumc+1 &&
                   jc > 1 && jc < xnumc && kc > 1 && kc < znumc+1
                    ΔRx_coarse[ic,jc,kc] = ΔRx_coarse[ic,jc,kc] / wtx[ic,jc,kc]
                else
                    ΔRx_coarse[ic,jc,kc] = 0
                end
            end
            # y-Stokes
            if ic < ynumc+1
                if wty[ic,jc,kc] != 0 && ic > 1 && ic < ynumc &&
                   jc > 1 && jc < xnumc+1 && kc > 1 && kc < znumc+1
                    ΔRy_coarse[ic,jc,kc] = ΔRy_coarse[ic,jc,kc] / wty[ic,jc,kc]
                else
                    ΔRy_coarse[ic,jc,kc] = 0
                end
            end
            # z-Stokes
            if kc < znumc+1
                if wtz[ic,jc,kc] != 0 && ic > 1 && ic < ynumc+1 &&
                   jc > 1 && jc < xnumc+1 && kc > 1 && kc < znumc
                    ΔRz_coarse[ic,jc,kc] = ΔRz_coarse[ic,jc,kc] / wtz[ic,jc,kc]
                else
                    ΔRz_coarse[ic,jc,kc] = 0
                end
            end
            # Continuity
            if ic < ynumc && jc < xnumc && kc < znumc
                if wtc[ic,jc,kc] != 0
                    ΔRc_coarse[ic,jc,kc] = ΔRc_coarse[ic,jc,kc] / wtc[ic,jc,kc] / etanc[ic,jc,kc]
                else
                    ΔRc_coarse[ic,jc,kc] = 0
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
    # Recomputing right parts (RX, RY, RZ, RC) for the coarser level (n+1).
    # Each kc-plane is independent (one write per coarse cell).
    parallel_kc = Threads.nthreads() > 1 && znumc >= RESTRICT_DIVIDE_K_PARALLEL_MIN_ZNUM
    if parallel_kc
        Threads.@threads for kc in 1:znumc+1
            _divide_restricted_residuals_on_coarser_kplane!(
                kc, xnumc, ynumc, znumc,
                ΔRx_coarse, ΔRy_coarse, ΔRz_coarse, ΔRc_coarse,
                wtx, wty, wtz, wtc, etanc,
            )
        end
    else
        for kc in 1:znumc+1
            _divide_restricted_residuals_on_coarser_kplane!(
                kc, xnumc, ynumc, znumc,
                ΔRx_coarse, ΔRy_coarse, ΔRz_coarse, ΔRc_coarse,
                wtx, wty, wtz, wtc, etanc,
            )
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
    
    @inbounds for j in 2:xnumf
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
    @inbounds for jc in 1:xnumc+1
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