module TimeStepCalculator

include("core/TimeStepUtils.jl")
include("core/TimeStepUtilsStepping.jl")
include("core/TimeStepExtension.jl")
include("core/TimeStepExtensionStepping.jl")

import EarthBox.ModelDataContainer: ModelData
import EarthBox.PrintFuncs: print_info
import EarthBox.ConversionFuncs: seconds_to_years, years_to_seconds
import EarthBox.BoundaryConditions: option_names, get_option_name
import .TimeStepExtension
import .TimeStepExtensionStepping

function update_time_step_parameters!(model::ModelData)::Nothing
    if use_velocity_stepping(model)
        update_time_step_parameters_with_stepping!(model)
    else
        update_time_step_parameters_no_stepping!(model)
    end
end

function update_time_step_parameters_no_stepping!(model::ModelData)::Nothing
    calculate_viscoelastic_time_step_no_stepping!(model)
    calculate_number_of_time_steps_no_stepping!(model)
    return nothing
end

function update_time_step_parameters_with_stepping!(
    model::ModelData,
)::Nothing
    
    _stepping_parameters_are_set = stepping_parameters_are_set(model)
    _velocity_parameters_are_set = velocity_parameters_are_set(model)
    using_symmetric_extension = check_for_symmetric_extension(model)

    if using_symmetric_extension && !_stepping_parameters_are_set
        if !timing_parameters_are_set(model)
            print_info("Required timing parameters: model_duration_myr, velocity_step_time", level=1)
            error("Timing parameters must be set for calculating time step parameters for velocity stepping.")
        end
        if !_velocity_parameters_are_set
            print_info("Required velocity parameters: full_velocity_extension, full_velocity_extension_step1", level=1)
            error("Velocity parameters must be set for calculating time step parameters for velocity stepping.")
        end
        TimeStepExtensionStepping.calculate_viscoelastic_time_steps_for_symmetric_extension_with_stepping(
            model,
            timestep_viscoelastic_limit_yrs = 50000.0,
            buffer_time_myr = 3.0
            )
    end
    return nothing
end

function timing_parameters_are_set(model::ModelData)::Bool
    model_duration_myr        = model.timestep.parameters.main_time_loop.model_duration_myr.value
    velocity_step_time        = model.bcs.parameters.velocity_step.velocity_step_time.value # Myr
    # Second step not included in this check as it is optional
    if isnan(model_duration_myr) || isnan(velocity_step_time)
        return false
    end
    return true
end

function velocity_parameters_are_set(model::ModelData)::Bool
    full_velocity_extension_m_s           = model.bcs.parameters.velocity.full_velocity_extension.value
    full_velocity_extension_m_s_step1     = model.bcs.parameters.velocity.full_velocity_extension_step1.value
    # Second step not included in this check as it is optional
    if isnan(full_velocity_extension_m_s) || isnan(full_velocity_extension_m_s_step1)
        return false
    end
    return true
end

function calculate_viscoelastic_time_step_no_stepping!(
    model::ModelData,
)::Nothing
    using_symmetric_extension = check_for_symmetric_extension(model)
    _is_time_step_set = is_time_step_set(model)
    if using_symmetric_extension && !_is_time_step_set
        (
            timestep_viscoelastic_s
        ) = TimeStepExtension.calculate_viscoelastic_time_step_for_symmetric_extension(model)
        print_info("Updated viscoelastic time step to account for symmetric extension", level=1)
        print_info("viscoelastic time step: $(seconds_to_years(timestep_viscoelastic_s)) yrs", level=2)
        main_time_loop = model.timestep.parameters.main_time_loop
        main_time_loop.timestep_viscoelastic.value = timestep_viscoelastic_s
    end
    return nothing
end

function is_time_step_set(model::ModelData)::Bool
    main_time_loop = model.timestep.parameters.main_time_loop
    timestep_viscoelastic = main_time_loop.timestep_viscoelastic.value
    if isnan(timestep_viscoelastic)
        return false
    end
    return true
end

function stepping_parameters_are_set(model::ModelData)::Bool
    main_time_loop = model.timestep.parameters.main_time_loop
    ntimestep_max = main_time_loop.ntimestep_max.value
    timestep_viscoelastic = main_time_loop.timestep_viscoelastic.value

    velocity_step_factor = model.bcs.parameters.velocity_step.velocity_step_factor.value
    are_not_loaded = isnan(ntimestep_max) || isnan(timestep_viscoelastic) || isnan(velocity_step_factor)
    if are_not_loaded
        return false
    end
    return true

end

function use_velocity_stepping(model::ModelData)::Bool
    iuse_velocity_step = model.bcs.parameters.velocity_step.iuse_velocity_step.value
    if iuse_velocity_step == 1
        return true
    end
    return false
end

function calculate_number_of_time_steps_no_stepping!(
    model::ModelData,
)::Nothing
    if check_for_symmetric_extension(model) && !is_number_of_time_steps_set(model)
        ntimestep_max = TimeStepExtension.calculate_number_of_time_steps_for_symmetric_extension(model)
        print_info("Updated number of time steps to account for symmetric extension", level=1)
        print_info("number of time steps: $(ntimestep_max)", level=2)
        main_time_loop = model.timestep.parameters.main_time_loop
        main_time_loop.ntimestep_max.value = ntimestep_max
    end
    return nothing
end

function is_number_of_time_steps_set(model::ModelData)::Bool
    main_time_loop = model.timestep.parameters.main_time_loop
    ntimestep_max = main_time_loop.ntimestep_max.value
    if ntimestep_max == -9999
        return false
    end
    return true
end

function update_transport_timestep_for_symmetric_extension!(
    model::ModelData, 
    number_of_transport_timesteps_per_model_timestep::Int
)::Nothing
    iuse_downhill_diffusion = model.obj_dict["iuse_downhill_diffusion"].value
    if iuse_downhill_diffusion == 1
        calculate_sediment_transport_timestep!(
            model, number_of_transport_timesteps_per_model_timestep)
    end
    return nothing
end

function calculate_sediment_transport_timestep!(
    model::ModelData,
    number_of_transport_timesteps_per_model_timestep::Int
)::Nothing
    using_sediment_transport = check_for_sediment_transport(model)
    if check_for_symmetric_extension(model) && !is_transport_time_step_set(model) && using_sediment_transport
        timestep_viscoelastic_s = model.timestep.parameters.main_time_loop.timestep_viscoelastic.value
        @assert !isnan(timestep_viscoelastic_s) "Timestep viscoelastic is NaN. Make sure you have set the viscoelastic time step."
        timestep_viscoelastic_yrs = seconds_to_years(timestep_viscoelastic_s)
        transport_timestep_yrs = (
            timestep_viscoelastic_yrs / Float64(number_of_transport_timesteps_per_model_timestep))
        transport_timestep_s = years_to_seconds(transport_timestep_yrs)
        model.topography.parameters.downhill_diffusion.transport_timestep.value = transport_timestep_s
        print_info("Updated transport time step to account for symmetric extension", level=1)
        print_info("Inputs:", level=2)
        print_info("timestep_viscoelastic_yrs: $(timestep_viscoelastic_yrs)", level=3)
        print_info("number_of_transport_timesteps_per_model_timestep: $(number_of_transport_timesteps_per_model_timestep)", level=3)
        print_info("Outputs:", level=2)
        print_info("transport time step: $(transport_timestep_yrs) yrs", level=2)
    end
    return nothing
end

function check_for_sediment_transport(model::ModelData)::Bool
    iuse_downhill_diffusion = model.topography.parameters.model_options.iuse_downhill_diffusion.value
    if iuse_downhill_diffusion == 1
        return true
    end
    return false
end

function is_transport_time_step_set(model::ModelData)::Bool
    transport_timestep = model.topography.parameters.downhill_diffusion.transport_timestep.value
    if isnan(transport_timestep)
        return false
    end
    return true
end

function check_for_symmetric_extension(model::ModelData)::Bool
    symmetric_extensional_models = [
        option_names.LithosphericExtensionMovingBoundaries,
        option_names.LithosphericExtensionFixedBoundaries,
        option_names.LithosphericExtensionDepthDependent,
        option_names.LithosphericExtensionInflowAndOutflowAlongSides,
    ]
    option_name = get_option_name(model)
    print_info("Checking for symmetric extension with option name: $(option_name)", level=1)
    return Symbol(option_name) in symmetric_extensional_models
end


end # module