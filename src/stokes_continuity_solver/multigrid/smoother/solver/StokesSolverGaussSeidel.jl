module StokesSolverGaussSeidel

include("solve/SolveStokes3dOpt.jl")
include("solve/SolveStokes2dOpt.jl")

import ...MultigridDataManager.MultigridStructures: RelaxationParameters
import ...LevelManager: LevelData, LevelData2d
import ...BoundaryConditions: set_boundary_conditions3d!
import ...BoundaryConditions: set_boundary_conditions2d!
import ...ArrayStats
import .SolveStokes3dOpt: solve_stokes_continuity_equations3d!
import .SolveStokes2dOpt: solve_stokes_continuity_equations2d!

const DEBUG = false

function solve_stokes_equations_gauss_seidel!(
    level_data::LevelData,
    smoothing_iterations::Vector{Int64},
    relaxation::RelaxationParameters
)::Nothing
    level_id = level_data.level_id
    nsmoothing_iterations = smoothing_iterations[level_id]
    for ismooth = 1:nsmoothing_iterations
        set_boundary_conditions3d!(level_data)
        solve_stokes_continuity_equations3d!(relaxation, level_data)
    end
    return nothing
end

function solve_stokes_equations_gauss_seidel!(
    level_data::LevelData2d,
    smoothing_iterations::Vector{Int64},
    relaxation::RelaxationParameters
)::Nothing
    level_id = level_data.level_id
    nsmoothing_iterations = smoothing_iterations[level_id]
    println(">>>> Solving stokes equations with Gauss-Seidel method for level $level_id with $nsmoothing_iterations iterations")
    dt_bc = 0.0
    dt_solve = 0.0
    for ismooth = 1:nsmoothing_iterations
        t1 = time()
        set_boundary_conditions2d!(level_data)
        t2 = time()
        dt_bc = dt_bc + t2-t1
        
        t1 = time()
        solve_stokes_continuity_equations2d!(relaxation, level_data)
        t2 = time()
        dt_solve = dt_solve + t2-t1
    end
    println(">>>> Total time taken to set boundary conditions: $(dt_bc) seconds")
    println(">>>> Total time taken to solve stokes equations: $(dt_solve) seconds")
    return nothing
end

end # module