"""
    MultigridSolver

Solution of Stokes and continuity equations in 3D with Multigrid based on V-cycle
by using external function Stokes_smoother3D(). Density distribution in the model 
corresponds to falling block test. Viscosity in the model is variable.

Staggered Grid for Multigrid

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
module MultigridSolver

import ..MultigridDataManager: MultigridData3d, MultigridData2d
import ..ViscosityScaling
import ..LevelManager
import ..LevelManager: LevelData, LevelData2d
import ..LevelInterpolation
import ..Smoother: update_rhs_parts_on_level1_using_global_level0_residuals
import ..MultigridVCycle
import ..ArrayStats

const DEBUG = true

function run_multigrid_solver(
    multigrid_data::Union{MultigridData3d, MultigridData2d}
)::Nothing
    println(">> Starting multigrid solver...")
    
    t1 = time()
    
    if multigrid_data.vcycle.use_multimulti
        update_rhs_parts_on_level1_using_global_level0_residuals(multigrid_data)
    end
    LevelInterpolation.interpolate_viscosity_to_coarser_levels!(multigrid_data)
    save_original_viscosity_on_all_layers!(multigrid_data)
    ViscosityScaling.rescale_viscosity!(multigrid_data)
    MultigridVCycle.execute_multigrid_vcycles!(multigrid_data)

    t2 = time()
    println(">> Finished multigrid solver in $(t2-t1) seconds")

    if DEBUG
        print_final_solution_stats(multigrid_data)
    end
    return nothing
end

function print_final_solution_stats(
    multigrid_data::Union{MultigridData3d, MultigridData2d}
)::Nothing
    if multigrid_data.vcycle.use_multimulti
        println(">> Final solution on level 0:")
        ArrayStats.print_solution_stats(multigrid_data.level0, "Final solution on level 0:")
    else
        println(">> Final solution on level 1:")
        ArrayStats.print_solution_stats(multigrid_data.level_vector[1], "Final solution on level 1:")
    end
    return nothing
end

function save_original_viscosity_on_all_layers!(
    multigrid_data::MultigridData3d
)::Nothing
    levelnum = length(multigrid_data.level_vector)
    for n = 1:levelnum
        LevelManager.save_original_viscosity!(multigrid_data.level_vector, n)
    end
    return nothing
end

function save_original_viscosity_on_all_layers!(
    multigrid_data::MultigridData2d
)::Nothing
    levelnum = length(multigrid_data.level_vector)
    for n = 1:levelnum
        LevelManager.save_original_viscosity!(multigrid_data.level_vector, n)
    end
    return nothing
end

end # module 