module SolveStokes3dOpt

import ....MultigridDataManager.MultigridStructures: RelaxationParameters
import ....LevelManager: LevelData, StencilArrays3d, current_stencil_arrays_3d
import ....Domain: on_vx_boundary3d_limits, on_vy_boundary3d_limits, on_vz_boundary3d_limits
import ....ArrayStats
import ...Residuals

# Parallelize over k only when there are enough k-planes; coarse grids pay thread overhead otherwise.
const GS_K_PARALLEL_MIN_ZNUM = 16

function _relax_gs_one_kplane!(
    color::Int,
    k::Int,
    Θ_stokes::Float64,
    Θ_continuity::Float64,
    xnum::Int,
    ynum::Int,
    znum::Int,
    A::StencilArrays3d,
)::Nothing
    ynum_vx = ynum + 1
    znum_vx = znum + 1
    xnum_vy = xnum + 1
    znum_vy = znum + 1
    xnum_vz = xnum + 1
    ynum_vz = ynum + 1
    # Red-black: for fixed (j,k), parity of i is determined by color. Stride i by 2 — half the visits,
    # no per-cell parity branch (same updates as scanning all i with a continue).
    @inbounds for j = 1:xnum+1
        i0 = (((1 + j + k) & 1) == color) ? 1 : 2
        for i in i0:2:ynum+1
            if j < xnum+1
                if !on_vx_boundary3d_limits(i, j, k, ynum_vx, znum_vx, xnum)
                    update_vx!(i, j, k, Θ_stokes, A)
                end
            end
            if i < ynum+1
                if !on_vy_boundary3d_limits(i, j, k, ynum, xnum_vy, znum_vy)
                    update_vy!(i, j, k, Θ_stokes, A)
                end
            end
            if k < znum+1
                if !on_vz_boundary3d_limits(i, j, k, ynum_vz, xnum_vz, znum)
                    update_vz!(i, j, k, Θ_stokes, A)
                end
            end
            if i < ynum && j < xnum && k < znum
                update_pressure!(i, j, k, Θ_continuity, A)
            end
        end
    end
    return nothing
end

function solve_stokes_continuity_equations3d!(
    relaxation::RelaxationParameters,
    level_data::LevelData
)::Nothing
    xnum = level_data.grid.parameters.geometry.xnum.value
    ynum = level_data.grid.parameters.geometry.ynum.value
    znum = level_data.grid.parameters.geometry.znum.value

    Θ_stokes = relaxation.relax_stokes
    Θ_continuity = relaxation.relax_continuity
    parallel_k = Threads.nthreads() > 1 && znum >= GS_K_PARALLEL_MIN_ZNUM
    A = current_stencil_arrays_3d(level_data)

    # Red-black (checkerboard) Gauss-Seidel: two sweeps per iteration.
    # Within each color, no two updated cells are direct neighbors,
    # so the sweep is safe for parallel execution.

    # (i + j + k) & 1 is the parity of the sum i + j + k: it is 0 when that sum 
    # is even, and 1 when it is odd. (Same idea as mod(i + j + k, 2), implemented with a bitwise AND.)
    
    # In Julia, & is bitwise AND on integers.
    # For two integers, each bit of the result is 1 only where both operands have a 1 in that position.
    # So (i + j + k) & 1 keeps only the least significant bit of the sum:
    # If the sum is even, that bit is 0.
    # If the sum is odd, that bit is 1.

    # The outer loop runs for color = 0:1 (lines 23–24), so there are two sweeps per relaxation pass:

    # color == 0: only grid points with even i + j + k are updated.
    # color == 1: only grid points with odd i + j + k are updated.
    # if (i + j + k) & 1 != color → continue means: “this (i,j,k) belongs to the other color; skip 
    # it on this sweep.”

    # So on a given color, the loop still visits every (i,j,k) in the nested ranges, but only 
    # executes the update_*! calls for cells on that color.
    for color = 0:1
        if parallel_k
            Threads.@threads for k = 1:znum+1
                _relax_gs_one_kplane!(
                    color, k, Θ_stokes, Θ_continuity, xnum, ynum, znum, A)
            end
        else
            for k = 1:znum+1
                _relax_gs_one_kplane!(
                    color, k, Θ_stokes, Θ_continuity, xnum, ynum, znum, A)
            end
        end
    end
    return nothing
end

# x-Stokes equation dSIGMAxx/dx+dSIGMAxy/dy+SIGMAxz/dz-dP/dx=RX
@inline function update_vx!(
    i::Int64,
    j::Int64,
    k::Int64,
    Θ_stokes::Float64,
    A::StencilArrays3d,
)::Nothing
    ΔR, Coef_vxC = Residuals.calculate_vx_residual(i, j, k, A)
    @inbounds A.vx[i,j,k] += ΔR/Coef_vxC*Θ_stokes
    return nothing
end

# y-Stokes equation dSIGMAyx/dx+dSIGMAyy/dy+SIGMAyz/dz-dP/dy=RY
@inline function update_vy!(
    i::Int64,
    j::Int64,
    k::Int64,
    Θ_stokes::Float64,
    A::StencilArrays3d,
)::Nothing
    ΔR, Coef_vyC = Residuals.calculate_vy_residual(i, j, k, A)
    @inbounds A.vy[i,j,k] += ΔR/Coef_vyC*Θ_stokes
    return nothing
end

# z-Stokes equation dSIGMAzx/dx+dSIGMAzy/dy+SIGMAzz/dz-dP/dz=RZ
@inline function update_vz!(
    i::Int64,
    j::Int64,
    k::Int64,
    Θ_stokes::Float64,
    A::StencilArrays3d,
)::Nothing
    ΔR, Coef_vzC = Residuals.calculate_vz_residual(i, j, k, A)
    @inbounds A.vz[i,j,k] += ΔR/Coef_vzC*Θ_stokes
    return nothing
end

@inline function update_pressure!(
    i::Int64,
    j::Int64,
    k::Int64,
    Θ_continuity::Float64,
    A::StencilArrays3d,
)::Nothing
    ΔR = Residuals.calculate_pressure_residual(i, j, k, A)
    @inbounds A.pr[i,j,k] += A.etan[i,j,k]*ΔR*Θ_continuity
    return nothing
end

end # module