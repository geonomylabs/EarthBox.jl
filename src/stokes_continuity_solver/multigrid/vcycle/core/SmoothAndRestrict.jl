module SmoothAndRestrict

import ...MultigridDataManager: MultigridData3d, MultigridData2d
import ...LevelManager: reset_solution_grids_to_zero!
import ...Restriction: restrict_stokes3d_residuals!
import ...Restriction: restrict_stokes2d_residuals!
import ...Smoother: stokes_continuity3d_viscous_smoother!
import ...Smoother: stokes_continuity2d_viscous_smoother!

function smooth_and_restrict!(
    multigrid_data::MultigridData3d,
)::Nothing
    pressure_bc = multigrid_data.pressure_bc
    smoothing_iterations = multigrid_data.smoothing_iterations
    relaxation = multigrid_data.relaxation
    level_vector = multigrid_data.level_vector

    levelnum = length(level_vector)
    for n = 1:levelnum
        if n > 1
            reset_solution_grids_to_zero!(level_vector[n])
        end
        (
            ΔRx_fine, ΔRy_fine, ΔRz_fine, ΔRc_fine
        ) = stokes_continuity3d_viscous_smoother!(
            pressure_bc, level_vector[n], smoothing_iterations, relaxation)
        if n < levelnum
            restrict_stokes3d_residuals!(
                n, level_vector, ΔRx_fine, ΔRy_fine, ΔRz_fine, ΔRc_fine)
        end
    end
    return nothing
end

function smooth_and_restrict!(
    multigrid_data::MultigridData2d,
)::Nothing
    pressure_bc = multigrid_data.pressure_bc
    smoothing_iterations = multigrid_data.smoothing_iterations
    relaxation = multigrid_data.relaxation
    level_vector = multigrid_data.level_vector

    levelnum = length(level_vector)
    for n = 1:levelnum
        if n > 1
            reset_solution_grids_to_zero!(level_vector[n])
        end
        (
            ΔRx_fine, ΔRy_fine, ΔRc_fine
        ) = stokes_continuity2d_viscous_smoother!(
            pressure_bc, level_vector[n], smoothing_iterations, relaxation)
        if n < levelnum
            restrict_stokes2d_residuals!(n, level_vector, ΔRx_fine, ΔRy_fine, ΔRc_fine)
        end
    end
    return nothing
end

end