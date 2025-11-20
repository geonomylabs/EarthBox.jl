"""
    SystemSolver

Module for solving system of equations for Stokes-continuity equations.
"""
module SystemSolver

import EarthBox.PrintFuncs: print_solution_vector_statistics
import SparseArrays: SparseMatrixCSC
import SparseArrays: sparse
import LinearAlgebra
import Printf: @printf
import EarthBox.PrintFuncs: @timeit_memit
import EarthBox.ModelDataContainer: ModelData
import EarthBox.Arrays: ArrayUtils
import EarthBox.ConfigurationManager.SolverConfig: SolverConfigState
import EarthBox.ParallelSolver: parallel_direct_solver
import ..GlobalIndices

"""
This function solves Stokes and Continuity equations on a 2D
fully staggered grid using a conservative finite difference approach to
discretize equations.
"""
function solve_system(
    model::ModelData,
    Li::Vector{Int64},
    Lj::Vector{Int64},
    Lv::Vector{Float64},
    solver_config::SolverConfigState
)::Tuple{Vector{Float64}, SparseMatrixCSC{Float64,Int64}}
    N = model.stokes_continuity.parameters.build.N.value
    Ls = build_sparse_crs_matrix(Li, Lj, Lv, N)
    use_mumps = solver_config.use_mumps
    if use_mumps
        R = copy(model.stokes_continuity.arrays.rhs.RHS.array)
        @timeit_memit "Finished parallel direct solver for Stokes-continuity equations" begin
            S = parallel_direct_solver(N, Li, Lj, Lv, R, solver_config)
        end
    else
        @timeit_memit "Finished serial direct solver for Stokes-continuity equations" begin
            S = stokes_solve_system_serial(model, Ls)
        end
    end
    print_solution_vector_statistics(S)
    return S, Ls
end

function build_sparse_crs_matrix(
    Li::Vector{Int64},
    Lj::Vector{Int64},
    Lv::Vector{Float64},
    N::Int64
)::SparseMatrixCSC{Float64,Int64}
    return sparse(Li, Lj, Lv, N, N)
end

function stokes_solve_system_serial(
    model::ModelData,
    Ls::SparseMatrixCSC{Float64,Int64}
)::Vector{Float64}
    rhs_arrays = model.stokes_continuity.arrays.rhs
    return Ls \ rhs_arrays.RHS.array
end

""" Initialize stokes-continuity parameters.

Steps:
1) Reset arrays solution arrays vx1, vy1 and pr1
2) Calculate pressure scaling coefficient
"""
function initialize_solver!(model::ModelData)::Nothing
    stokes_arrays = model.stokes_continuity.arrays
    ArrayUtils.setzeros!(stokes_arrays.staggered_grid_velocity.vx1)
    ArrayUtils.setzeros!(stokes_arrays.staggered_grid_velocity.vy1)
    ArrayUtils.setzeros!(stokes_arrays.pressure.pr1)

    stokes_params = model.stokes_continuity.parameters
    stokes_params.build.pscale.value = pressure_scaling_coeff(model)
    return nothing
end

function pressure_scaling_coeff(model::ModelData)::Float64
    xstpavg = model.grids.parameters.geometry.xstpavg.value
    ystpavg = model.grids.parameters.geometry.ystpavg.value
    eta = model.stokes_continuity.arrays.viscosity.etan0.array[1, 1]
    return 2.0 * eta / (xstpavg + ystpavg)
end

function stokes_solve_system!(
    model::ModelData,
    Ls::SparseMatrixCSC{Float64,Int64}
)::Nothing
    rhs_arrays = model.stokes_continuity.arrays.rhs
    S = Ls \ rhs_arrays.RHS.array
    copy_solution_array_to_model_data_structure!(model, S)
    return nothing
end

function stokes_calc_min_max_solu(model::ModelData)::Nothing
    solu_arrays = model.stokes_continuity.arrays.stokes_solution
    Smin = minimum(solu_arrays.soluv1.array)
    Smax = maximum(solu_arrays.soluv1.array)
    @printf("min and max solution vec %f %f\n", Smin, Smax)
    return nothing
end

""" Process Stokes-Continuity solution vector and large-matrix.

Process Stokes solution vector including loading velocity solutions
and boundary conditions into staggered grids arrays vx1 and vy1 and
load pressure solution into pressure grid array pr1.
"""
function process_stokes_solution!(
    model::ModelData,
    S::Vector{Float64}
)::Nothing
    @timeit_memit "Finished processing Stokes-continuity solution" begin
        copy_solution_array_to_model_data_structure!(model, S)
        reload_solu_vel!(model)
        apply_vel_bc!(model)
    end
    return nothing
end

""" Copy solution array to model data structure.

Updated Arrays:
    model.stokes_continuity.arrays.stokes_solution
    ----------------------------------------------
    soluv1: Vector{Float64}
        Stokes-continuity solution vector with vx, vy and pr solutions
"""
function copy_solution_array_to_model_data_structure!(
    model::ModelData,
    S::Vector{Float64}
)::Nothing
    model.stokes_continuity.arrays.stokes_solution.soluv1.array = copy(S)
    return nothing
end

""" Load stokes-continuity solution to velocity and pressure arrays.

The solution vector is transferred to staggered velocity grid arrays
and pressure grid array.

Updated Arrays:
    model.stokes_continuity.arrays.staggered_grid_velocity
    ------------------------------------------------------
    vx1.array: Matrix{Float64}
        X-component of velocity (m/s) defined on staggered vx grid.

    vy1.array: Matrix{Float64}
        Y-component of velocity (m/s) defined on staggered vy grid.

    model.stokes_continuity.arrays.pressure
    ---------------------------------------
    pr1.array: Matrix{Float64}
        Pressure (Pa) on pressure grid.
"""
function reload_solu_vel!(model::ModelData)::Nothing
    ArrayUtils.setzeros!(model.stokes_continuity.arrays.staggered_grid_velocity.vy1)
    ArrayUtils.setzeros!(model.stokes_continuity.arrays.staggered_grid_velocity.vx1)
    ArrayUtils.setzeros!(model.stokes_continuity.arrays.pressure.pr1)

    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value
    pscale = model.stokes_continuity.parameters.build.pscale.value
    S = model.stokes_continuity.arrays.stokes_solution.soluv1.array
    vx1 = model.stokes_continuity.arrays.staggered_grid_velocity.vx1.array
    vy1 = model.stokes_continuity.arrays.staggered_grid_velocity.vy1.array
    pr1 = model.stokes_continuity.arrays.pressure.pr1.array

    for j in 1:(xnum-1)
        for i in 1:(ynum-1)
            # Indexes for P, vx and vy
            cell_index = GlobalIndices.get_global_basic_cell_index(i, j, ynum)
            ivx = GlobalIndices.get_global_ivx_unknown_index(cell_index)
            ivy = GlobalIndices.get_global_ivy_unknown_index(ivx)
            ipr = GlobalIndices.get_global_ipr_unknown_index(ivx)
            # Reload Vx
            if j < xnum - 1
                vx1[i+1, j+1] = S[ivx]
            end
            # Reload Vy
            if i < ynum-1
                vy1[i+1, j+1] = S[ivy]
            end
            # Reload P
            pr1[i, j] = S[ipr]*pscale
        end
    end
    return nothing
end

"""
    apply_vel_bc(model::ModelData)

Apply velocity boundary conditions to staggered velocity grids.

Updated Arrays:
    model.stokes_continuity.arrays.staggered_grid_velocity
    ------------------------------------------------------
    vx1.array: Matrix{Float64}
        X-component of velocity (m/s) defined on staggered vx grid.

    vy1.array: Matrix{Float64}
        Y-component of velocity (m/s) defined on staggered vy grid.
"""
function apply_vel_bc!(model::ModelData)::Nothing
    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value
    bleftx = model.bcs.arrays.vel_comp.bleftx.array
    blefty = model.bcs.arrays.vel_comp.blefty.array
    brightx = model.bcs.arrays.vel_comp.brightx.array
    brighty = model.bcs.arrays.vel_comp.brighty.array
    btopx = model.bcs.arrays.vel_comp.btopx.array
    btopy = model.bcs.arrays.vel_comp.btopy.array
    bbottomx = model.bcs.arrays.vel_comp.bbottomx.array
    bbottomy = model.bcs.arrays.vel_comp.bbottomy.array
    vx1 = model.stokes_continuity.arrays.staggered_grid_velocity.vx1.array
    vy1 = model.stokes_continuity.arrays.staggered_grid_velocity.vy1.array

    # Top and bottom boundaries
    for j in 1:xnum
        vx1[1, j] = btopx[j, 1] + btopx[j, 2]*vx1[2, j]
        vx1[ynum+1, j] = bbottomx[j, 1] + bbottomx[j, 2]*vx1[ynum, j]
    end
    for j in 1:(xnum+1)
        vy1[1, j] = btopy[j, 1] + btopy[j, 2]*vy1[2, j]
        vy1[ynum, j] = bbottomy[j, 1] + bbottomy[j, 2]*vy1[ynum-1, j]
    end

    # Left and right boundaries
    for i in 1:(ynum+1)
        vx1[i, 1] = bleftx[i, 1] + bleftx[i, 2]*vx1[i, 2]
        vx1[i, xnum] = brightx[i, 1] + brightx[i, 2]*vx1[i, xnum-1]
    end
    for i in 1:ynum
        vy1[i, 1] = blefty[i, 1] + blefty[i, 2]*vy1[i, 2]
        vy1[i, xnum+1] = brighty[i, 1] + brighty[i, 2]*vy1[i, xnum]
    end
    return nothing
end

"""
    backup_stokes_arrays(model::ModelData)

Backup Stokes solution arrays.

Updated Arrays:
    model.stokes_continuity.arrays.stokes_solution
    ----------------------------------------------
    soluv1_old.array: Vector{Float64}
        Copy of array soluv1, Stokes-continuity solution vector with vx, vy
        and pr solutions.

    model.stokes_continuity.arrays.staggered_grid_velocity
    ------------------------------------------------------
    vx1_old.array: Matrix{Float64}
        Copy of array vx1, x-component of velocity (m/s) defined on staggered
        vx grid.

    vy1_old.array: Matrix{Float64}
        Copy of array vy1, y-component of velocity (m/s) defined on staggered
        vy grid.

    model.stokes_continuity.arrays.pressure
    ---------------------------------------
    pr1_old.array: Matrix{Float64}
        Copy of array pr1, pressure (Pa) on pressure grid.

    model.stokes_continuity.arrays.velocity_solution
    ------------------------------------------------
    vxy_old.array: Vector{Float64}
       Copy of vxy array, solution array with vx and vy solutions
       (no pressure).
"""
function backup_stokes_arrays!(model::ModelData)::Nothing
    solu_arrays = model.stokes_continuity.arrays.stokes_solution
    sgvel_arrays = model.stokes_continuity.arrays.staggered_grid_velocity
    pr_arrays = model.stokes_continuity.arrays.pressure
    velsolu_arrays = model.stokes_continuity.arrays.velocity_solution

    solu_arrays.soluv1_old.array = copy(solu_arrays.soluv1.array)
    sgvel_arrays.vx1_old.array = copy(sgvel_arrays.vx1.array)
    sgvel_arrays.vy1_old.array = copy(sgvel_arrays.vy1.array)
    pr_arrays.pr1_old.array = copy(pr_arrays.pr1.array)
    velsolu_arrays.vxy_old.array = copy(velsolu_arrays.vxy.array)
    return nothing
end

end # module 