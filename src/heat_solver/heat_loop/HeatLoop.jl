module HeatLoop

import Printf: @sprintf
import EarthBox.PrintFuncs: @timeit_memit, print_info
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ConfigurationManager.SolverConfig: SolverConfigState
import ..SystemSolver: discretize_and_solve_conductive_heat_equation
import ..SystemSolver: HeatResiduals

function run_conductive_heat_solver_loop!(
    model::ModelData,
    solver_config::SolverConfigState
)::Nothing
    """ Solve conductive temperature equation using adaptive time stepping.

    A new temperature solution array `tk2` on the basic grid is obtained by 
    solving the heat conduction equation using adaptive time stepping if the 
    maximum temperature change from array `dtk1` exceeds a specified limit 
    (max_temp_change). If adaptive time stepping is activated, a loop is used to
    solve the heat equation until the total time (timestep_sum) is equal to the
    main model time step (timestep_start). The old temperature solution array 
    `tk0`, which defines the initial conditions for each iteration, is 
    initialized by copying the temperature array `tk1` that was previously 
    interpolated from markers. The new temperature array `tk2` is copied to 
    `tk0` to define the old solution during heat loop iterations. Note that 
    when the loop terminates `tk0` is equal to `tk2` since `tk2` is copied
    to tk0 at the end of the loop sequence.

    This function also updates the temperature residuals array `rest` during
    each iteration of the heat loop.

    # Arguments
    - `model::ModelData`: The model data container object used to store 
        parameters and arrays.
    - `solver_config::SolverConfig`: Configuration object containing the solver 
        settings.

    # Updated Arrays
    - `model.heat_equation.arrays.temperature.tk0`: Old temperature (K) solution 
        on basic grid
    - `model.heat_equation.arrays.temperature.tk2`: Updated conductive 
        temperature solution on basic grid
    - `model.heat_equation.arrays.temperature.rest`: Temperature residuals
    """
    sec_per_myr = model.conversion.parameters.sec_per_Myr.value
    max_temp_change_limit = model.heat_equation.parameters.temp_change_limit.max_temp_change.value

    timestep_start = model.timestep.parameters.main_time_loop.timestep.value
    timestep_loop = timestep_start
    set_thermal_loop_time_step(model, timestep_loop)

    timestep_sum = 0.0
    set_total_thermal_loop_time(model, timestep_sum)

    icount = 0
    while timestep_sum < timestep_start
        print_info(model, icount, timestep_sum, timestep_start, sec_per_myr)
        discretize_and_solve_conductive_heat_equation(model, solver_config)
        HeatResiduals.compute_residuals_temp!(model)
        update_temperature_change_array(model)
        dtkmax = calculate_maximum_temperature_change(model)
        #println("    => dtkmax: $dtkmax max_temp_change_limit: $max_temp_change_limit")
        if dtkmax > max_temp_change_limit
            timestep_loop = timestep_loop * max_temp_change_limit / dtkmax
            #println("    => timestep_loop: $timestep_loop")
            set_thermal_loop_time_step(model, timestep_loop)
            print_heat_loop_info(model)
            discretize_and_solve_conductive_heat_equation(model, solver_config)
            HeatResiduals.compute_residuals_temp!(model)
        end
        timestep_sum = timestep_sum + timestep_loop
        set_total_thermal_loop_time(model, timestep_sum)
        timestep_loop = min(timestep_loop, timestep_start - timestep_sum)
        set_thermal_loop_time_step(model, timestep_loop)
        copy_new_temperature_to_old(model)
        icount = icount + 1
    end
    return nothing
end

function print_info(
    model::ModelData,
    icount::Int,
    timestep_sum::Float64,
    timestep_start::Float64,
    sec_per_myr::Float64
)::Nothing
    print_iteration_info(icount, timestep_sum, timestep_start, sec_per_myr)
    print_heat_loop_info(model)
    return nothing
end

function print_iteration_info(
    icount::Int,
    timestep_sum::Float64,
    timestep_start::Float64,
    sec_per_myr::Float64
)::Nothing
    print_info(
        "Heat loop iteration $icount :" *
        "timestep_sum (Myr): $(timestep_sum/sec_per_myr) :" *
        "timestep_start (Myr): $(timestep_start/sec_per_myr)",
        level=2
    )
    return nothing
end

function update_temperature_change_array(model::ModelData)::Nothing
    temperature = model.heat_equation.arrays.temperature
    temperature.dtk1.array .= temperature.tk2.array .- temperature.tk0.array
    return nothing
end

@inline function calculate_maximum_temperature_change(model::ModelData)::Float64
    return maximum(abs.(model.heat_equation.arrays.temperature.dtk1.array))
end

@inline function set_thermal_loop_time_step(model::ModelData, timestep_loop::Float64)::Nothing
    # This is the time step used to solve the transient heat conduction equation.
    model.timestep.parameters.thermal_loop.timestep_heat.value = timestep_loop
    return nothing
end

function print_heat_loop_info(model::ModelData)::Nothing
    sec_per_myr = model.conversion.parameters.sec_per_Myr.value

    timestep = model.timestep.parameters.main_time_loop.timestep.value
    timestep_heat = model.timestep.parameters.thermal_loop.timestep_heat.value
    timestep_sum = model.timestep.parameters.thermal_loop.timestep_sum.value

    ts_start = timestep/sec_per_myr
    ts_current = timestep_heat/sec_per_myr
    ts_post_solution = timestep_sum/sec_per_myr + timestep_heat/sec_per_myr

    print_info(
        "Target Time Step (Myr): $(@sprintf("%.8f", ts_start))" *
        ": Current Time Step (Myr): $(@sprintf("%.8f", ts_current))" *
        ": Integrated Time Step (Myr): $(@sprintf("%.8f", ts_post_solution))",
        level=3
    )
    return nothing
end

@inline function set_total_thermal_loop_time(model::ModelData, timestep_sum::Float64)::Nothing
    model.timestep.parameters.thermal_loop.timestep_sum.value = timestep_sum
    return nothing
end

function copy_new_temperature_to_old(model::ModelData)::Nothing
    model.heat_equation.arrays.temperature.tk0.array .= 
        model.heat_equation.arrays.temperature.tk2.array
    return nothing
end

end # module 