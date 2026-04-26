module StokesResiduals

include("utils/ResidualStructs.jl")
include("utils/Continuity.jl")
include("utils/XStokes.jl")
include("utils/YStokes.jl")

import LinearAlgebra
import LinearAlgebra: mul!, norm
import EarthBox.ModelDataContainer: ModelData
import EarthBox.Arrays: ArrayUtils
import SparseArrays: SparseMatrixCSC
import .XStokes
import .YStokes
import .Continuity
import .ResidualStructs: StokesResidualsInput

function stokes_calc_residual_L2norm(model::ModelData, Ls::Any)::Float64
    work = model.stokes_continuity.arrays.residuals.resnl_work
    mul!(work, Ls, model.stokes_continuity.arrays.stokes_solution.soluv1.array)
    work .-= model.stokes_continuity.arrays.rhs.RHS.array
    return norm(work)
end

""" Calculate non-linear residual.

# Updated array from group model.stokes_continuity.arrays.residuals:
- `resnl.array::Vector{Float64}` ((xnum-1)*(ynum-1)*3):
    - Non-linear stokes-continuity system residual.
"""
function stokes_calc_nonlinear_system_residual!(
    model::ModelData,
    Ls::SparseMatrixCSC{Float64,Int64}
)::Nothing
    work = model.stokes_continuity.arrays.residuals.resnl_work
    resnl = model.stokes_continuity.arrays.residuals.resnl.array
    mul!(work, Ls, model.stokes_continuity.arrays.stokes_solution.soluv1_old.array)
    work .-= model.stokes_continuity.arrays.rhs.RHS.array
    map!(abs, resnl, work)
    return nothing
end

""" Allocation-free variant of `stokes_calc_nonlinear_system_residual!`.

Computes the nonlinear residual `|A * x_old - R|` directly from the
preallocated COO triplets in
`model.stokes_continuity.parameters.build.system_vectors` instead of
materializing a `SparseMatrixCSC` first.

The result is mathematically equal to the legacy variant — both compute
`sum(Lv[k] * x_old[Lj[k]] for k where Li[k] == i)` for each row `i`, with
duplicate `(i, j)` triplets summed naturally. Floating-point results
agree to machine-epsilon (relative ~1e-16) but are not bit-identical:
`sparse()` reorders the triplets into CSC layout and `mul!` then sums
column-major, while this loop sums in original triplet order, so the FP
summation order differs. The Picard outer-loop convergence behavior is
unaffected by this rounding difference.

# Updated array from group model.stokes_continuity.arrays.residuals:
- `resnl.array::Vector{Float64}` ((xnum-1)*(ynum-1)*3):
    - Non-linear stokes-continuity system residual.
"""
function stokes_calc_nonlinear_system_residual_optimized!(
    model::ModelData
)::Nothing
    work = model.stokes_continuity.arrays.residuals.resnl_work
    resnl = model.stokes_continuity.arrays.residuals.resnl.array
    sv = model.stokes_continuity.parameters.build.system_vectors
    Li = sv.Li_out
    Lj = sv.Lj_out
    Lv = sv.Lv_out
    x_old = model.stokes_continuity.arrays.stokes_solution.soluv1_old.array
    R = model.stokes_continuity.arrays.rhs.RHS.array
    nnz = length(Li)
    fill!(work, 0.0)
    @inbounds for k in 1:nnz
        work[Li[k]] += Lv[k] * x_old[Lj[k]]
    end
    work .-= R
    map!(abs, resnl, work)
    return nothing
end

""" Calculate stokes-continuity residuals.

# Updated arrays from group model.stokes_continuity.arrays.residuals:
- `resx1::Matrix{Float64}` ((ynum+1), (xnum)): 
    - Residual for x-Stokes equation on staggered vx grid.
- `resy1::Matrix{Float64}` ((ynum), (xnum+1)): 
    - Residual for y-Stokes equation on staggered vy grid.
- `resc1::Matrix{Float64}` ((ynum-1), (xnum-1)): 
    - Residual for continuity equation on pressure grid.
"""
function compute_stokes_residuals!(model::ModelData)::Nothing
    ArrayUtils.setzeros!(model.stokes_continuity.arrays.residuals.resx1)
    ArrayUtils.setzeros!(model.stokes_continuity.arrays.residuals.resy1)
    ArrayUtils.setzeros!(model.stokes_continuity.arrays.residuals.resc1)

    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value

    residuals_input = StokesResidualsInput(
        xnum,
        ynum,
        model.grids.arrays.basic.xstp_b.array,
        model.grids.arrays.basic.ystp_b.array,
        model.grids.arrays.staggered_vy.xstp_vy.array,
        model.grids.arrays.staggered_vx.ystp_vx.array,
        model.bcs.arrays.internal.bintern_zone.array,
        model.stokes_continuity.arrays.viscosity.etas0.array,
        model.stokes_continuity.arrays.viscosity.etan0.array,
        model.stokes_continuity.arrays.rhs.RX1.array,
        model.stokes_continuity.arrays.rhs.RY1.array,
        model.stokes_continuity.arrays.rhs.RC1.array,
        model.stokes_continuity.arrays.staggered_grid_velocity.vx1.array,
        model.stokes_continuity.arrays.staggered_grid_velocity.vy1.array,
        model.stokes_continuity.arrays.pressure.pr1.array
    )

    Threads.@threads for j in 1:(xnum+1)
        for i in 1:(ynum+1)
            XStokes.calculate_residual!(
                residuals_input, i, j,
                model.stokes_continuity.arrays.residuals.resx1.array
            )
            YStokes.calculate_residual!(
                residuals_input, i, j,
                model.stokes_continuity.arrays.residuals.resy1.array
            )
            Continuity.calculate_residual!(
                residuals_input, i, j,
                model.stokes_continuity.arrays.residuals.resc1.array
            )
        end
    end
    return nothing
end

""" Calculate non-linear Stokes-continuity residuals.

# Updated Arrays
- `model.stokes_continuity.arrays.residuals.resnlx`: Residual for x-Stokes 
  equation on staggered vx grid.
- `model.stokes_continuity.arrays.residuals.resnly`: Residual for y-Stokes 
  equation on staggered vy grid.
- `model.stokes_continuity.arrays.residuals.resnlc`: Residual for continuity 
  equation on pressure grid.
"""
function compute_stokes_nonlinear_residuals!(model::ModelData)::Nothing
    ArrayUtils.setzeros!(model.stokes_continuity.arrays.residuals.resnlx)
    ArrayUtils.setzeros!(model.stokes_continuity.arrays.residuals.resnly)
    ArrayUtils.setzeros!(model.stokes_continuity.arrays.residuals.resnlc)

    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value

    residuals_input = StokesResidualsInput(
        xnum,
        ynum,
        model.grids.arrays.basic.xstp_b.array,
        model.grids.arrays.basic.ystp_b.array,
        model.grids.arrays.staggered_vy.xstp_vy.array,
        model.grids.arrays.staggered_vx.ystp_vx.array,
        model.bcs.arrays.internal.bintern_zone.array,
        model.stokes_continuity.arrays.viscosity.etas0.array,
        model.stokes_continuity.arrays.viscosity.etan0.array,
        model.stokes_continuity.arrays.rhs.RX1.array,
        model.stokes_continuity.arrays.rhs.RY1.array,
        model.stokes_continuity.arrays.rhs.RC1.array,
        model.stokes_continuity.arrays.staggered_grid_velocity.vx1_old.array,
        model.stokes_continuity.arrays.staggered_grid_velocity.vy1_old.array,
        model.stokes_continuity.arrays.pressure.pr1_old.array
    )

    Threads.@threads for j in 1:(xnum+1)
        for i in 1:(ynum+1)
            XStokes.calculate_residual!(
                residuals_input, i, j,
                model.stokes_continuity.arrays.residuals.resnlx.array
            )
            YStokes.calculate_residual!(
                residuals_input, i, j,
                model.stokes_continuity.arrays.residuals.resnly.array
            )
            Continuity.calculate_residual!(
                residuals_input, i, j,
                model.stokes_continuity.arrays.residuals.resnlc.array
            )
        end
    end
    return nothing
end

end # module 