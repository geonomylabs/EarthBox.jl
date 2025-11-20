""" Iteratively solve non-linear Stokes-continuity equations.

This module uses a node-based global loop over Stokes-continuity equations.
The Picard method is used to solve the non-linear system.
"""
module PlasticityLoopNodes

import EarthBox.PrintFuncs: print_info
import EarthBox.ModelDataContainer: ModelData
import EarthBox: StokesContinuitySolver
import EarthBox.Rheology.PlasticFailure.PlasticFailureNodes: 
    update_nodes_for_plastic_yielding
import EarthBox.Interpolation: PlasticStrainRateInterpolation
import EarthBox: Interpolation
import EarthBox.LoopInputStruct: LoopInput
import ...Initialize
import ...LoopInformation
import ...LoopTermination
import ...UpdateLoopParameters
import ...Convergence

"""
    picard_loop(model::ModelData, loop_input::LoopInput)

Picard loop with node-based plasticity.

This iteratively solves the non-linear Stokes-continuity equations using
a Picard iteration method using node-based plasticity. The loop terminates 
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
)::Nothing
    stokes_solution_norms = model.stokes_continuity.parameters.solution_norms

    tol_global = model.stokes_continuity.parameters.picard.tolerance_picard.value
    nglobal = model.stokes_continuity.parameters.picard.nglobal.value

    Initialize.initialize_stokes_global_loop_parameters!(model)

    _iglobal = 0
    for _iglobal in 1:nglobal
        if _iglobal > 1
            Interpolation.interpolate_basic_viscoplastic_viscosity_to_pressure_grid!(model)
        end
        StokesContinuitySolver.solve_viscoelastic_stokes_continuity_equations!(
            model, loop_input.solver_config)
        StokesContinuitySolver.stokes_solution_processing!(model, inside_flags)

        (
            etavp_new,
            plastic_strain_rate,
            plastic_yield_new,
            nnode_yield,
            global_yield_error
        ) = update_nodes_for_plastic_yielding(
            model, loop_input.no_yielding_in_mobile_wall,
            loop_input.no_yielding_in_plate_extension
        )

        if nnode_yield > 0
            update_viscoplastic_viscosity_on_basic_grid_for_plastic_failure!(
                model, etavp_new)
            update_plastic_yield_flags_on_basic_grid!(
                model, plastic_yield_new)
            update_plastic_strain_rate_invariant_on_basic_and_pressure_grids!(
                model, plastic_strain_rate)
        end
        LoopInformation.process_global_loop_info_and_update_dict!(
            model, _iglobal, nglobal, global_yield_error)
        UpdateLoopParameters.update!(model)
        convergence_criterion = stokes_solution_norms.dvxy_rel_L2.value

        stop_plastic_loop = LoopTermination.check_stop_global_loop_nodes(
            nnode_yield, _iglobal, nglobal,
            global_yield_error, convergence_criterion, tol_global
        )
        if stop_plastic_loop
            break
        end
    end
    Convergence.save_iteration_info_to_file(
        loop_input.run_paths["output_dir"], model.global_iter_dict)
    return nothing
end

function print_empty_space()
    print_info("", level=1)
end

function calc_global_yield_error(
    global_yield_error_sum::Float64,
    nnode_yield::Int
)::Float64
    global_yield_error = sqrt(global_yield_error_sum / Float64(nnode_yield))
    return global_yield_error
end

function update_basic_grid_plasticity!(
    model::ModelData,
    etavp_new::Matrix{Float64},
    plastic_yield_new::Matrix{Float64}
)
    model.stokes_continuity.arrays.viscosity.etas1.array = copy(etavp_new)
    model.stokes_continuity.arrays.plastic_def.plastic_yield.array = copy(plastic_yield_new)
end

function update_plastic_yield_flags_on_basic_grid!(
    model::ModelData,
    plastic_yield_new::Matrix{Float64}
)
    model.stokes_continuity.arrays.plastic_def.plastic_yield.array = copy(plastic_yield_new)
end

function update_viscoplastic_viscosity_on_basic_grid_for_plastic_failure!(
    model::ModelData,
    etavp_new::Matrix{Float64}
)
    model.stokes_continuity.arrays.viscosity.etas1.array = copy(etavp_new)
end

function update_plastic_strain_rate_invariant_on_basic_and_pressure_grids!(
    model::ModelData,
    plastic_strain_rate::Matrix{Float64}
)
    strain_rate_and_spin = model.stokes_continuity.arrays.strain_rate_and_spin
    strain_rate_and_spin.eii_plastic_basic.array = copy(plastic_strain_rate)
    PlasticStrainRateInterpolation.interpolate_at_pressure_nodes!(model)
end

function reset_grid_plasticity!(
    model::ModelData,
    etas1_initial::Matrix{Float64},
    plastic_yield_initial::Matrix{Float64}
)
    model.stokes_continuity.arrays.viscosity.etas1.array = copy(etas1_initial)
    model.stokes_continuity.arrays.plastic_def.plastic_yield.array = copy(plastic_yield_initial)
end

function decrease_timestep!(model::ModelData)
    dtkoef = 2.0
    timestep = model.timestep.parameters.main_time_loop.timestep.value
    timestep = timestep / dtkoef
    model.timestep.parameters.main_time_loop.timestep.value = timestep
end

end # module 