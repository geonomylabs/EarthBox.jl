module TimeStepExtensionStepping

import EarthBox.PrintFuncs: print_info
import EarthBox.ModelDataContainer: ModelData
import EarthBox: OptionTools
import EarthBox.StaggeredGrid: check_for_ttype_refinement
import EarthBox.ConversionFuncs: years_to_seconds, seconds_to_years, m_s_to_cm_yr
import ..TimeStepUtilsStepping: get_time_step_parameters_for_velocity_stepping
import ..TimeStepUtils: calculate_minimum_cell_spacing

""" 
     calculate_viscoelastic_time_steps_for_symmetric_extension_with_stepping(
        model::ModelData, 
        timestep_viscoelastic_limit_yrs::Float64, 
        buffer_time_myr::Float64
    )::Nothing

Time step and stepping parameters are calculated using stepping timing and
full extension velocity at each step.

# Arguments
- `model::ModelData`
    - Model data container containing the model parameters and arrays
- `timestep_viscoelastic_limit_yrs::Float64`
    - Time step viscoelastic limit in years
- `buffer_time_myr::Float64`
    - Buffer time in millions of years added to model duration when calculating 
       the maximum number of time steps

# Returns
- `nothing`

"""
function calculate_viscoelastic_time_steps_for_symmetric_extension_with_stepping(
    model::ModelData;
    timestep_viscoelastic_limit_yrs::Float64,
    buffer_time_myr::Float64
)::Nothing
    model_duration_myr        = model.timestep.parameters.main_time_loop.model_duration_myr.value
    velocity_step_time        = model.bcs.parameters.velocity_step.velocity_step_time.value # Myr
    velocity_second_step_time = model.bcs.parameters.velocity_step.velocity_second_step_time.value # Myr

    full_velocity_extension_m_s           = model.bcs.parameters.velocity.full_velocity_extension.value
    full_velocity_extension_m_s_step1     = model.bcs.parameters.velocity.full_velocity_extension_step1.value
    full_velocity_extension_m_s_step2     = model.bcs.parameters.velocity.full_velocity_extension_step2.value
    characteristic_max_velocity_m_s       = full_velocity_extension_m_s / 2.0
    characteristic_max_velocity_m_s_step1 = full_velocity_extension_m_s_step1 / 2.0

    if !isnan(full_velocity_extension_m_s_step2)
        characteristic_max_velocity_m_s_step2 = full_velocity_extension_m_s_step2 / 2.0
    else
        characteristic_max_velocity_m_s_step2 = NaN
    end

    marker_cell_displ_max   = model.markers.parameters.advection.marker_cell_displ_max.value
    minimum_cell_spacing_m  = calculate_minimum_cell_spacing(model)

    velocity_step_params = get_time_step_parameters_for_velocity_stepping(
        model_duration_myr,
        velocity_step_time, # Myr 
        velocity_second_step_time, # Myr
        characteristic_max_velocity_m_s,
        characteristic_max_velocity_m_s_step1,
        characteristic_max_velocity_m_s_step2,
        minimum_cell_spacing_m,
        marker_cell_displ_max,
        timestep_viscoelastic_limit_yrs,
        buffer_time_myr
        )
        
    ntimestep_max                     = velocity_step_params["ntimestep_max"]
    velocity_step_factor              = velocity_step_params["velocity_step_factor"]
    velocity_second_step_factor       = velocity_step_params["velocity_second_step_factor"]
    timestep_adjustment_factor        = velocity_step_params["timestep_adjustment_factor"]
    timestep_second_adjustment_factor = velocity_step_params["timestep_second_adjustment_factor"]
    timestep_viscoelastic             = velocity_step_params["timestep_viscoelastic"]
    timestep_viscoelastic_step1       = velocity_step_params["timestep_viscoelastic_step1"]
    timestep_viscoelastic_step2       = velocity_step_params["timestep_viscoelastic_step2"]

    main_time_loop = model.timestep.parameters.main_time_loop
    main_time_loop.ntimestep_max.value               = ntimestep_max
    main_time_loop.timestep_viscoelastic.value       = timestep_viscoelastic
    main_time_loop.timestep_viscoelastic_step1.value = timestep_viscoelastic_step1
    main_time_loop.timestep_viscoelastic_step2.value = timestep_viscoelastic_step2

    velocity_step = model.bcs.parameters.velocity_step
    velocity_step.velocity_step_factor.value              = velocity_step_factor
    velocity_step.velocity_second_step_factor.value       = velocity_second_step_factor
    velocity_step.timestep_adjustment_factor.value        = timestep_adjustment_factor
    velocity_step.timestep_second_adjustment_factor.value = timestep_second_adjustment_factor

    print_info("Calculated time step and velocity stepping parameters for symmetric extension", level=1)
    print_info("Inputs:", level=2)
    print_info("model_duration: $(model_duration_myr) (Myr)", level=3)
    print_info("velocity_step_time: $(velocity_step_time) (Myr)", level=3)
    print_info("velocity_second_step_time: $(velocity_second_step_time) (Myr)", level=3)
    print_info("characteristic_max_velocity: $(m_s_to_cm_yr(characteristic_max_velocity_m_s)) (cm/yr)", level=3)
    print_info("characteristic_max_velocity_step1: $(m_s_to_cm_yr(characteristic_max_velocity_m_s_step1)) (cm/yr)", level=3)
    print_info("characteristic_max_velocity_step2: $(m_s_to_cm_yr(characteristic_max_velocity_m_s_step2)) (cm/yr)", level=3)
    print_info("minimum_cell_spacing_m: $(minimum_cell_spacing_m) (m)", level=3)
    print_info("marker_cell_displ_max: $(marker_cell_displ_max) (m)", level=3)
    print_info("timestep_viscoelastic_limit: $(timestep_viscoelastic_limit_yrs) (yrs)", level=3)
    print_info("Outputs Used to Update Model):", level=2)
    print_info("ntimestep_max: $(ntimestep_max)", level=3)
    print_info("timestep_viscoelastic: $(seconds_to_years(timestep_viscoelastic)) (yrs)", level=3)
    print_info("timestep_viscoelastic_step1: $(seconds_to_years(timestep_viscoelastic_step1)) (yrs)", level=3)
    print_info("timestep_viscoelastic_step2: $(seconds_to_years(timestep_viscoelastic_step2)) (yrs)", level=3)
    print_info("velocity_step_factor: $(velocity_step_factor)", level=3)
    print_info("velocity_second_step_factor: $(velocity_second_step_factor)", level=3)
    print_info("timestep_adjustment_factor: $(timestep_adjustment_factor)", level=3)
    print_info("timestep_second_adjustment_factor: $(timestep_second_adjustment_factor)", level=3)
    return nothing
end

end # module