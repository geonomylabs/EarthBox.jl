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

"""When `ENV["EARTHBOX_MG_DEBUG"] != "1"`, skip verbose final solution statistics."""
mg_debug_stats_enabled() = (get(ENV, "EARTHBOX_MG_DEBUG", "1") == "1")


"""Run the multigrid solver.

Here's a summary of all the optimizations Cursor(Opus 4.6) implemented:

Phase 1: Reduce Per-Iteration Cost

@inbounds on all hot-loop kernels - Added to SolveStokes3dOpt.jl, SolveStokes2dOpt.jl, Residuals.jl 
(all 4 residual calculation functions + the compute loop), Restriction.jl (all 4 loop nests), 
Prolongation.jl (both 2D/3D loops), TriLinearInterpolation.jl, and BiLinearInterpolation.jl. Also 
cached scalar_fine[i,j,k] to a local val in interpolation functions to avoid redundant lookups.

Fused viscosity scaling - Replaced exp.(log.(...)) with @. fused broadcasting using (x/min_o)^eta_coeff, 
eliminating all temporary array allocations in scale_viscosity_for_all_levels!.

Pre-allocated residual arrays - Added res_vx_buf, res_vy_buf, res_vz_buf, res_pr_buf to LevelData. 
compute_residuals! now fills them in-place with fill! + reuse instead of allocating zeros(...) every call.

Pre-allocated restriction workspace - Restriction now writes directly into level_vector[n+1].RX.array etc. 
instead of allocating new arrays. Weight arrays restrict_wtx/wty/wtz/wtc are pre-allocated in LevelData. 
The etan_resc temp is pre-allocated in etan_resc_buf and filled with @. fusion.

Pre-allocated prolongation corrections - prolongate_stokes3d_solution now uses prolong_dvx/dvy/dvz/dpr 
buffers from the finer LevelData instead of allocating via grid_array3D.

In-place prolongation relaxation - Replaced dvx * relax_velocity temporaries with @. fused operations: @. 
finer.vx.array += dvx * relax_velocity.

Phase 2: Reduce Total V-Cycle Count

Multimulti convergence (3D) uses max_global_residual only (aligned with 2D multimulti). Principle
residuals are still accumulated each V-cycle for diagnostics/plots but are not used for stopping,
since they can show transient sawtooth drops that do not reflect sustained convergence.

Smoothing cap on all coarse levels - set_smoothing_iterations now applies max_smoothing_on_coarsest 
as a cap to ALL levels (not just the last), preventing over-smoothing on intermediate coarse grids 
(e.g., 80 iterations on a 7^3 grid is now capped).

Red-black Gauss-Seidel ordering - Split the single threaded GS sweep into two passes: "red" cells 
(i+j+k) % 2 == 0 then "black" cells (i+j+k) % 2 == 1. Each pass is safe for parallel execution since 
no two cells of the same color are direct neighbors, eliminating data races from the previous approach.

Profiling (one-off)

Set `ENV["EARTHBOX_MG_TIMING"]="1"` and optionally `ENV["EARTHBOX_MG_TIMING_DETAIL"]="1"` for a breakdown
of smoother vs restriction vs prolongation inside the restrict/prolong legs (see MultigridVCycle).

For a sampling profile in the REPL after `using Profile` and importing your driver (e.g. StokesSinker):

    Profile.clear()
    @profile StokesSinker.run_stokes_sinker(make_plots=false, use_multimulti=true, model_type=:ThreeDimensional)
    Profile.print(maxdepth=30)

Use a shorter grid or fewer `nvcycles` if the run is too long. For flame graphs, use ProfileView.jl or
record with `Profile.Allocs` if investigating allocations.

3D restriction (weight divide on the coarse grid) and prolongation parallelize over k-planes when
`Threads.nthreads() > 1` and the relevant z extent is at least 16 (same threshold style as GS/residuals).

Optional algorithmic next steps (not implemented here; benchmark before/after)

- Direct or sparse solve on the coarsest staggered level instead of many GS sweeps.
- W-cycle or F-cycle (more coarse work per outer cycle; may reduce outer iteration count).
- `restot` in ViscosityScaling.jl (non-multimulti ramp) uses `resx` four times; using all four principle
  components changes continuation timing and is a tuning decision, not a pure micro-optimization.

"""
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

    if mg_debug_stats_enabled()
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