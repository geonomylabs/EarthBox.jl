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
    for ismooth = 1:nsmoothing_iterations
        set_boundary_conditions2d!(level_data)
        solve_stokes_continuity_equations2d!(relaxation, level_data)
    end
    return nothing
end

end # module