""" Iteratively solve non-linear Stokes-continuity equations.

This module uses a marker-based global loop over Stokes-continuity equations.
The Picard method is used to solve the non-linear system.
"""
module PlasticityLoopMarkers

import EarthBox.PrintFuncs: print_info
import EarthBox.ModelDataContainer: ModelData
import EarthBox: StokesContinuitySolver
import EarthBox.Interpolation: BackupTransportArrays
import EarthBox.Interpolation: MarkersToStokesArrays
import EarthBox.Rheology.PlasticFailure: MarkerPlasticity
import EarthBox: Interpolation
import EarthBox.LoopInputStruct: LoopInput
import ...Initialize
import ...LoopInformation
import ...LoopTermination
import ...UpdateLoopParameters
import ...Convergence

"""
    picard_loop(model::ModelData, loop_input::LoopInput)

Picard loop with marker-based plasticity.

This iteratively solves the non-linear Stokes-continuity equations using
a Picard iteration method using marker-based plasticity. The loop terminates 
when the convergence criterion is met or the maximum number of iterations is
reached.

# Arguments
- `model::ModelData`: Model data container
- `loop_input::LoopInput`: Loop input struct
"""
function picard_loop(
    model::ModelData,
    loop_input::LoopInput,
    inside_flags::Vector{Int8}
)
    stokes_solution_norms = model.stokes_continuity.parameters.solution_norms

    tol_global = model.stokes_continuity.parameters.picard.tolerance_picard.value
    nglobal = model.stokes_continuity.parameters.picard.nglobal.value
    no_yielding_in_mobile_wall = loop_input.no_yielding_in_mobile_wall

    Initialize.initialize_stokes_global_loop_parameters!(model)

    global_yield_error = 1e38  # Not used for marker-based global loop
    _iglobal = 0
    for _iglobal in 1:nglobal
        if _iglobal > 1
            BackupTransportArrays.backup_and_clear_viscosity_transport_arrays!(model)
            MarkerPlasticity.update_plastic_yielding!(model, inside_flags, no_yielding_in_mobile_wall)
            MarkersToStokesArrays.interpolate_viscosity_to_transport_arrays!(model, inside_flags)
            Interpolation.interpolate_basic_viscoplastic_viscosity_to_pressure_grid!(model)
        end
        StokesContinuitySolver.solve_viscoelastic_stokes_continuity_equations!(
            model, loop_input.solver_config)
        StokesContinuitySolver.stokes_solution_processing!(model, inside_flags)
        LoopInformation.process_global_loop_info_and_update_dict!(
            model, _iglobal, nglobal, global_yield_error)
        UpdateLoopParameters.update!(model)
        convergence_criterion = stokes_solution_norms.dvxy_rel_L2.value
        stop_plastic_loop = LoopTermination.check_stop_global_loop_markers(
            _iglobal, nglobal, convergence_criterion, tol_global)
        if stop_plastic_loop
            break
        end
    end
    Convergence.save_iteration_info_to_file(
        loop_input.run_paths["output_dir"], model.global_iter_dict)
end

function print_empty_space()
    print_info("", level=1)
end

end # module 