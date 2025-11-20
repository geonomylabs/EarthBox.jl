module SystemSolver

include("residuals/HeatResiduals.jl")
include("build_system/HeatBuildManager.jl")

import EarthBox.PrintFuncs: print_solution_vector_statistics
import SparseArrays: SparseMatrixCSC
import SparseArrays: sparse
import Printf
import EarthBox.PrintFuncs: @timeit_memit
import EarthBox.BuildSysTools: SystemVectors
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ConfigurationManager.SolverConfig: SolverConfigState
import EarthBox.ParallelSolver: parallel_direct_solver
import .HeatBuildManager

""" Build and solve system of equations for the 2D conductive heat equation.

This function manages the building and solving of the system of equations
for the 2D heat conduction equation. The 2D conductive heat equation is
discretized on the irregularly spaced basic grid using a conservative 
finite-difference approximation. Right-hand side terms are defined using 
array RT defined on the basic grid. Thermal Boundary condition are defined 
using arrays bleftt, brightt, btopt and bbott. The initial temperature is
defined by array tk0. The final temperature solution is copied to array tk2.

Updated Array Objects
-------------------
tk2: Matrix{Float64}
    Updated temperature solution on basic grid.
"""
function discretize_and_solve_conductive_heat_equation(
    model::ModelData,
    solver_config::SolverConfigState
)::Nothing
    @timeit_memit "Finished building the system of equations for heat equation" begin
        system_vectors = HeatBuildManager.build_system_of_equations(model)
    end
    N = model.heat_equation.parameters.build.N.value
    Ls = build_sparse_crs_matrix(system_vectors, N)
    use_mumps = solver_config.use_mumps
    
    if use_mumps
        R = copy(model.heat_equation.arrays.rhs.RHSheat.array)
        @timeit_memit "Finished parallel direct solver for heat equation" begin
            S = parallel_direct_solver(
                N, system_vectors.Li, system_vectors.Lj, system_vectors.Lv, 
                R, solver_config
            )
        end
    else
        @timeit_memit "Finished serial direct solver for heat equation" begin
            S = heat_solve_system_serial(model, Ls)
        end
    end
    #print_solution_vector_statistics(S)
    reload_solu_temp(model, S)
    return nothing
end

function build_sparse_crs_matrix(
    system_vectors::SystemVectors,
    N::Int64
)::SparseMatrixCSC{Float64,Int64}
    return sparse(system_vectors.Li, system_vectors.Lj, system_vectors.Lv, N, N)
end

function heat_solve_system_serial(
    model::ModelData,
    Ls::SparseMatrixCSC{Float64,Int64}
)::Vector{Float64}
    rhs = model.heat_equation.arrays.rhs
    return Ls \ rhs.RHSheat.array
end

function reload_solu_temp(
    model::ModelData,
    S::Vector{Float64}
)::Nothing
    grid_geom = model.grids.parameters.geometry
    temperature = model.heat_equation.arrays.temperature
    reload_solu_temp_loop!(
        grid_geom.xnum.value, grid_geom.ynum.value, S, temperature.tk2.array)
    return nothing
end

function reload_solu_temp_loop!(
    xnum::Int64,
    ynum::Int64,
    S::Vector{Float64},
    tk2::Matrix{Float64}
)::Nothing
    for j in 1:xnum
        for i in 1:ynum
            # Index for T
            itk = (j-1)*ynum + i
            # Reload T
            tk2[i, j] = S[itk]
        end
    end
    return nothing
end

end # module 