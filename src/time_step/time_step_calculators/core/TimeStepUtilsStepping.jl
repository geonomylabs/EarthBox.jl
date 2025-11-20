module TimeStepUtilsStepping

import EarthBox.PrintFuncs: print_info
import EarthBox.ConversionFuncs: years_to_seconds, seconds_to_years
import ..TimeStepUtils: calculate_viscoelastic_time_step

""" Calculate time step parameters for velocity stepping.

# Arguments
- `case_inputs::CaseType`: Dictionary of case inputs with units converted to 
  standard units.
- `print_output::Bool`: Whether to print output information (default: true)

# Returns
- `velocity_step_params::Dict`: Dictionary of velocity step parameters with 
  the following keys:
  - ntimestep_max: Maximum number of time steps
  - velocity_step_factor: Velocity step factor
  - velocity_second_step_factor: Velocity second step factor
  - timestep_adjustment_factor: Time step adjustment factor
  - timestep_second_adjustment_factor: Time step second adjustment factor
  - timestep_viscoelastic_yrs: Initial viscoelastic time step in years
  - timestep_viscoelastic_yrs_step1: Viscoelastic time step in years for step 1
  - timestep_viscoelastic_yrs_step2: Viscoelastic time step in years for step 2
"""
function get_time_step_parameters_for_velocity_stepping(
    model_duration_myr::Float64, # Myr
    velocity_step_time::Float64, # Myr
    velocity_second_step_time::Float64, # Myr
    characteristic_max_velocity_m_s::Float64,
    characteristic_max_velocity_m_s_step1::Float64,
    characteristic_max_velocity_m_s_step2::Float64,
    minimum_cell_spacing_m::Float64,
    marker_cell_displ_max::Float64,
    timestep_viscoelastic_limit_yrs::Float64,
    buffer_time_myr::Float64
)::Dict{String, Union{Float64, Int}}
    
    (
        timestep_viscoelastic_yrs, 
        timestep_viscoelastic_yrs_step1,
        timestep_viscoelastic_yrs_step2
    ) = calculate_viscoelastic_time_steps(
        characteristic_max_velocity_m_s,
        characteristic_max_velocity_m_s_step1,
        characteristic_max_velocity_m_s_step2,
        minimum_cell_spacing_m,
        marker_cell_displ_max,
        timestep_viscoelastic_limit_yrs
        )

    velocity_step_factor = characteristic_max_velocity_m_s_step1 / characteristic_max_velocity_m_s
    if !isnan(characteristic_max_velocity_m_s_step2)
        velocity_second_step_factor = characteristic_max_velocity_m_s_step2 / characteristic_max_velocity_m_s_step1
    else
        velocity_second_step_factor = NaN
    end

    timestep_adjustment_factor = timestep_viscoelastic_yrs_step1 / timestep_viscoelastic_yrs
    if !isnan(timestep_viscoelastic_yrs_step2)
        timestep_second_adjustment_factor = timestep_viscoelastic_yrs_step2 / timestep_viscoelastic_yrs_step1
    else
        timestep_second_adjustment_factor = NaN
    end

    ntimestep_max = calculate_total_number_of_time_steps_for_velocity_step_model(
        timestep_viscoelastic_yrs,
        timestep_viscoelastic_yrs_step1,
        timestep_viscoelastic_yrs_step2,
        velocity_step_time,
        velocity_second_step_time,
        model_duration_myr,
        buffer_time_myr
        )

    velocity_step_params = Dict{String, Union{Float64, Int}}(
        "ntimestep_max"                     => ntimestep_max,
        "velocity_step_factor"              => velocity_step_factor,
        "velocity_second_step_factor"       => velocity_second_step_factor,
        "timestep_adjustment_factor"        => timestep_adjustment_factor,
        "timestep_second_adjustment_factor" => isnan(timestep_second_adjustment_factor) ? NaN : timestep_second_adjustment_factor,
        "timestep_viscoelastic"             => years_to_seconds(timestep_viscoelastic_yrs),
        "timestep_viscoelastic_step1"       => years_to_seconds(timestep_viscoelastic_yrs_step1),
        "timestep_viscoelastic_step2"       => isnan(timestep_viscoelastic_yrs_step2) ? NaN : years_to_seconds(timestep_viscoelastic_yrs_step2)
    )

    return velocity_step_params
end

function calculate_viscoelastic_time_steps(
    characteristic_max_velocity_m_s::Float64,
    characteristic_max_velocity_m_s_step1::Float64,
    characteristic_max_velocity_m_s_step2::Float64,
    minimum_cell_spacing_m::Float64,
    marker_cell_displ_max::Float64,
    timestep_viscoelastic_limit_yrs::Float64
)::Tuple{Float64, Float64, Float64}
    timestep_viscoelastic_yrs = calculate_viscoelastic_time_step(
        characteristic_max_velocity_m_s,
        minimum_cell_spacing_m,
        marker_cell_displ_max,
        timestep_viscoelastic_limit_yrs
    )
    timestep_viscoelastic_yrs_step1 = calculate_viscoelastic_time_step(
        characteristic_max_velocity_m_s_step1,
        minimum_cell_spacing_m,
        marker_cell_displ_max,
        timestep_viscoelastic_limit_yrs
    )
    timestep_viscoelastic_yrs_step2 = calculate_viscoelastic_time_step(
        characteristic_max_velocity_m_s_step2,
        minimum_cell_spacing_m,
        marker_cell_displ_max,
        timestep_viscoelastic_limit_yrs
    )
    return (
        timestep_viscoelastic_yrs, 
        timestep_viscoelastic_yrs_step1,
        timestep_viscoelastic_yrs_step2
        )
end

function calculate_total_number_of_time_steps_for_velocity_step_model(
    timestep_viscoelastic_yrs_initial::Float64,
    timestep_viscoelastic_yrs_step1::Float64,
    timestep_viscoelastic_yrs_step2::Float64,
    velocity_step_time_myr::Float64,
    velocity_second_step_time_myr::Float64,
    model_duration_myr::Float64,
    buffer_time_myr::Float64
)::Int
    number_of_timesteps_initial = calculate_number_of_timesteps(
        timestep_viscoelastic_yrs_initial, velocity_step_time_myr)
    if isnan(velocity_second_step_time_myr)
        delta_time_myr = model_duration_myr - velocity_step_time_myr
        buffer_time_myr_first_step = buffer_time_myr
    else
        delta_time_myr = velocity_second_step_time_myr - velocity_step_time_myr
        # Do not use buffer time for first step if second step is defined
        buffer_time_myr_first_step = 0.0
    end
    number_of_timesteps_step1 = calculate_number_of_timesteps(
        timestep_viscoelastic_yrs_step1, 
        delta_time_myr, 
        buffer_time_myr = buffer_time_myr_first_step
        )
    if has_second_step_inputs(timestep_viscoelastic_yrs_step2, velocity_second_step_time_myr)
        number_of_timesteps_step2 = calculate_number_of_timesteps(
            timestep_viscoelastic_yrs_step2, 
            model_duration_myr - velocity_second_step_time_myr;
            buffer_time_myr = buffer_time_myr
        )
    else
        number_of_timesteps_step2 = 0.0
    end
    total_number_of_time_steps = (
        number_of_timesteps_initial
        + number_of_timesteps_step1
        + number_of_timesteps_step2
        )
    return total_number_of_time_steps
end

function has_second_step_inputs(
    timestep_viscoelastic_yrs_step2::Union{Float64, Nothing},
    velocity_second_step_time_myr::Union{Float64, Nothing},
)::Bool
    has_inputs = true
    if isnan(timestep_viscoelastic_yrs_step2) || isnothing(timestep_viscoelastic_yrs_step2)
        has_inputs = false
    end
    if isnan(velocity_second_step_time_myr) || isnothing(velocity_second_step_time_myr)
        has_inputs = false
    end
    return has_inputs
end

function calculate_number_of_timesteps(
    timestep_viscoelastic_yrs::Float64,
    duration_myr::Float64;
    buffer_time_myr::Float64=0.0
)::Int
    ntimestep_max = floor(Int, duration_myr * 1e6 / timestep_viscoelastic_yrs)
    if buffer_time_myr > 0.0
        ntimestep_max += floor(Int, buffer_time_myr * 1e6 / timestep_viscoelastic_yrs)
    end
    return ntimestep_max
end

""" Calculate number of output time steps for velocity stepping.

If iend is defined then the number of time steps is calculated as the
difference between iend and istart plus one taking into account a one 
starting point.

# Arguments
- `case_inputs::CaseType`: Dictionary of case inputs with units converted to 
  standard units.
- `istart::Int`: Index of the first time step (default: 1)
- `iend::Union{Int, Nothing}`: Index of the last time step (default: nothing)

# Returns
- `number_of_output_time_steps::Float64`: Number of output time steps
"""
function calculate_number_of_output_time_steps_for_velocity_stepping(
    model_duration_myr::Float64,
    velocity_step_time::Float64, # Myr
    velocity_second_step_time::Float64, # Myr
    timestep_viscoelastic_yrs::Float64,
    timestep_viscoelastic_yrs_step1::Float64,
    timestep_viscoelastic_yrs_step2::Float64,
    timestep_out_yrs::Float64,
    istart::Int=1,
    iend::Union{Int, Nothing}=nothing
)::Float64
    if iend === nothing
        number_of_output_time_step1 = calculate_incremental_output_time_steps(
            velocity_step_time * 1e6,
            timestep_viscoelastic_yrs,
            timestep_out_yrs
            )
        number_of_output_time_step2 = calculate_incremental_output_time_steps(
            (velocity_second_step_time - velocity_step_time) * 1e6,
            timestep_viscoelastic_yrs_step1, 
            timestep_out_yrs
            )
        number_of_output_time_step3 = calculate_incremental_output_time_steps(
            (model_duration_myr - velocity_second_step_time) * 1e6,
            timestep_viscoelastic_yrs_step2,
            timestep_out_yrs
            )
        ntimestep_max = (
            number_of_output_time_step1
            + number_of_output_time_step2
            + number_of_output_time_step3
        )
    else
        ntimestep_max = (iend - istart) + 1
    end
    return ntimestep_max
end

function calculate_incremental_output_time_steps(
    model_duration_yrs::Float64,
    timestep_viscoelastic_yrs::Float64,
    timestep_out_yrs::Float64
)::Int
    number_of_output_time_steps = floor(Int,
        model_duration_yrs / timestep_viscoelastic_yrs /
        floor(Int, timestep_out_yrs / timestep_viscoelastic_yrs))
    return number_of_output_time_steps
end

function calculate_velocity_step_adjustment_factors(
    characteristic_max_velocity::Float64,
    characteristic_max_velocity_step1::Float64,
    characteristic_max_velocity_step2::Float64
)::Tuple{Float64, Float64}
    velocity_step1_factor = characteristic_max_velocity_step1 / characteristic_max_velocity
    velocity_step2_factor = characteristic_max_velocity_step2 / characteristic_max_velocity_step1
    return velocity_step1_factor, velocity_step2_factor
end

end