module TimeStepExtension

import EarthBox.PrintFuncs: print_info
import EarthBox.ModelDataContainer: ModelData
import EarthBox: OptionTools
import EarthBox.StaggeredGrid: check_for_ttype_refinement
import EarthBox.ConversionFuncs: years_to_seconds, seconds_to_years
import ..TimeStepUtils: calculate_viscoelastic_time_step
import ..TimeStepUtils: calculate_number_of_timesteps
import ..TimeStepUtils: calculate_minimum_cell_spacing

function calculate_viscoelastic_time_step_for_symmetric_extension(
    model::ModelData,
    timestep_viscoelastic_limit_yrs::Float64 = 50000.0
)::Float64
    full_velocity_extension_m_s     = model.bcs.parameters.velocity.full_velocity_extension.value
    marker_cell_displ_max           = model.markers.parameters.advection.marker_cell_displ_max.value
    minimum_cell_spacing_m          = calculate_minimum_cell_spacing(model)
    characteristic_max_velocity_m_s = full_velocity_extension_m_s / 2.0

    timestep_viscoelastic_yrs = calculate_viscoelastic_time_step(
        characteristic_max_velocity_m_s, 
        minimum_cell_spacing_m, 
        marker_cell_displ_max, 
        timestep_viscoelastic_limit_yrs
        )
    print_info("Inputs:", level=2)
    print_info("full_velocity_extension_m_s: $(full_velocity_extension_m_s)", level=3)
    print_info("marker_cell_displ_max: $(marker_cell_displ_max)", level=3)
    print_info("minimum_cell_spacing_m: $(minimum_cell_spacing_m)", level=3)
    print_info("characteristic_max_velocity_m_s: $(characteristic_max_velocity_m_s)", level=3)
    print_info("timestep_viscoelastic_limit_yrs: $(timestep_viscoelastic_limit_yrs)", level=3)
    print_info("Outputs:", level=2)
    print_info("timestep_viscoelastic_yrs: $(timestep_viscoelastic_yrs)", level=3)

    return years_to_seconds(timestep_viscoelastic_yrs)
end

function calculate_number_of_time_steps_for_symmetric_extension(
    model::ModelData,
    buffer_time_myr::Float64 = 3.0
)::Int
    model_duration_myr = model.timestep.parameters.main_time_loop.model_duration_myr.value
    timestep_viscoelastic_s = model.timestep.parameters.main_time_loop.timestep_viscoelastic.value
    timestep_viscoelastic_yrs = seconds_to_years(timestep_viscoelastic_s)
    ntimestep_max = calculate_number_of_timesteps(
        timestep_viscoelastic_yrs, model_duration_myr, buffer_time_myr
        )
    return ntimestep_max
end

end # module