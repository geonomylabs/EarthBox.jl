module VelocityStep

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: load_parameters!
import EarthBox.ModelDataContainer: set_parameter!
import EarthBox.InitializationTools: update_use_flags
import EarthBox.ConversionFuncs: seconds_to_years
import EarthBox.ConversionFuncs: meters_per_seconds_to_cm_per_yr

const PDATA = get_eb_parameters()

struct ValidInputNames
    iuse_velocity_step::Symbol
    velocity_step_factor::Symbol
    timestep_adjustment_factor::Symbol
    timestep_second_adjustment_factor::Symbol
    velocity_step_time::Symbol
    velocity_second_step_time::Symbol
    velocity_second_step_factor::Symbol
end

"""
    initialize!(model::ModelData; kwargs...)::Nothing

Initialize velocity stepping model. Activating velocity stepping requires setting the 
`iuse_velocity_step` parameter to 1. At the velocity step time the lateral
velocity boundary conditions are multiplied by the `velocity_step_factor` 
parameter. This is accompanied by a multiplying the viscoelastic time step (`timestep_viscoelastic`) by 
the `timestep_adjustment_factor` used to maintain stability during an increase in 
velocity. There are two available velocity steps. The second velocity step is
ignored if associated parameters are not provided.

For certain cases involving symmetric extensional boundary conditions with two velocity stepping events, 
the time step parameters, velocity stepping factors and time step adjustment factors are calculated
automatically. See [`EarthBox.run_time_steps`](@ref EarthBox.run_time_steps) for more details.
For these cases, ensure that velocity stepping factors and time step adjustment factors are not
overridden by providing values for the `velocity_step_factor`, `velocity_second_step_factor`,
`timestep_adjustment_factor` and `timestep_second_adjustment_factor` parameters. The integer flag
`iuse_velocity_step` should be still be set to 1 to activate velocity stepping and the first and second 
velocity step times (`velocity_step_time` and `velocity_second_step_time`) should be provided.

# Arguments
- `model::`[`ModelData`](@ref ModelData): Model data object containing model parameters and arrays.

# Optional Keyword Arguments
- `iuse_velocity_step::Int64`: 
    - $(PDATA.iuse_velocity_step.description)
- `velocity_step_factor::Float64`: 
    - $(PDATA.velocity_step_factor.description)
- `velocity_second_step_factor::Float64`: 
    - $(PDATA.velocity_second_step_factor.description)
- `timestep_adjustment_factor::Float64`: 
    - $(PDATA.timestep_adjustment_factor.description)
- `timestep_second_adjustment_factor::Float64`: 
    - $(PDATA.timestep_second_adjustment_factor.description)
- `velocity_step_time::Float64`: 
    - $(PDATA.velocity_step_time.description)
- `velocity_second_step_time::Float64`: 
    - $(PDATA.velocity_second_step_time.description)

# Returns
- `Nothing`

# Example: Single velocity step
``` julia
using EarthBox
eb = EarthBoxState(
    xnum = 100, ynum = 100, 
    xsize = 1000000.0, ysize = 1000000.0, 
    dx_marker = 10000.0, dy_marker = 10000.0
)
model = eb.model_manager.model
BoundaryConditions.VelocityStep.initialize!(
    model, iuse_velocity_step = 1, velocity_step_factor = 2.0, 
    timestep_adjustment_factor = 0.5, velocity_step_time = 10.0
)
```

"""
function initialize!(model::ModelData; kwargs...)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    return nothing
end

function reset_bc_for_veloc_step!(model::ModelData)::Nothing
    iuse = model.bcs.parameters.velocity_step.iuse_velocity_step.value
    if iuse == 1 && at_step_time(model)
        println(">> Resetting extensional/contractional velocity for velocity step.")
        update_extension_velocity!(model)
        update_contraction_velocity!(model)
        update_timestep_viscoelastic!(model)
        update_nskip!(model)
        update_ivelocity_step_counter!(model)
    end
    return nothing
end

function at_step_time(model::ModelData)::Bool
    timesum = model.timestep.parameters.main_time_loop.timesum.value
    timesum_myr = seconds_to_years(timesum)/1e6

    ivelocity_step_counter = model.bcs.parameters.velocity_step.ivelocity_step_counter.value
    velocity_first_step_time, velocity_second_step_time = get_velocity_step_times(model)

    if ivelocity_step_counter in (0, 1)
        velocity_step_time = ivelocity_step_counter == 0 ? velocity_first_step_time : velocity_second_step_time
        if velocity_step_time > 0.0 # Disregard default values
            return timesum_myr >= velocity_step_time
        end
    end
    return false
end

function print_step_info(
    ivelocity_step_counter::Int,
    velocity_first_step_time::Float64,
    velocity_second_step_time::Float64,
    velocity_step_time::Float64,
    check::Bool
)::Nothing
    println(">> ivelocity_step_counter: ", ivelocity_step_counter)
    println(">> First velocity step time (Myr): ", velocity_first_step_time)
    println(">> Second velocity step time (Myr): ", velocity_second_step_time)
    println(">> Current velocity step time (Myr): ", velocity_step_time)
    println(">> At velocity step time: ", check)
    return nothing
end

function get_velocity_step_times(model::ModelData)::Tuple{Float64, Float64}
    velocity_first_step_time = model.bcs.parameters.velocity_step.velocity_step_time.value
    velocity_second_step_time = model.bcs.parameters.velocity_step.velocity_second_step_time.value
    return (velocity_first_step_time, velocity_second_step_time)
end

function update_extension_velocity!(model::ModelData)::Nothing
    velocity_step_factor = get_step_velocity_factor(model)
    full_velocity_extension = model.bcs.parameters.velocity.full_velocity_extension.value
    full_velocity_extension *= velocity_step_factor

    println(
        ">> New full_velocity_extension after step (cm/yr): ",
        meters_per_seconds_to_cm_per_yr(full_velocity_extension)
    )
    model.bcs.parameters.velocity.full_velocity_extension.value = full_velocity_extension
    return nothing
end

function update_contraction_velocity!(model::ModelData)::Nothing
    velocity_step_factor = get_step_velocity_factor(model)
    full_velocity_contraction = model.bcs.parameters.velocity.full_velocity_contraction.value
    full_velocity_contraction *= velocity_step_factor

    println(
        ">> New full_velocity_contraction after step (cm/yr): ",
        meters_per_seconds_to_cm_per_yr(full_velocity_contraction)
    )
    model.bcs.parameters.velocity.full_velocity_contraction.value = full_velocity_contraction
    return nothing
end

function get_step_velocity_factor(model::ModelData)::Float64
    ivelocity_step_counter = model.bcs.parameters.velocity_step.ivelocity_step_counter.value
    velocity_step = model.bcs.parameters.velocity_step

    if ivelocity_step_counter == 0
        return velocity_step.velocity_step_factor.value
    elseif ivelocity_step_counter == 1
        return velocity_step.velocity_second_step_factor.value
    else
        return 1.0
    end
end

function update_timestep_viscoelastic!(model::ModelData)::Nothing
    main_time_loop = model.timestep.parameters.main_time_loop
    timestep_viscoelastic = main_time_loop.timestep_viscoelastic.value
    timestep_adjustment_factor = get_timestep_adjustment_factor(model)

    if abs(timestep_adjustment_factor) != 1.0
        timestep_viscoelastic *= timestep_adjustment_factor
        main_time_loop.timestep_viscoelastic.value = timestep_viscoelastic
        println(
            ">> New time step (yr) after velocity step: ",
            seconds_to_years(timestep_viscoelastic)
        )
    end
    return nothing
end

function update_nskip!(model::ModelData)::Nothing
    timestep_adjustment_factor = get_timestep_adjustment_factor(model)

    if abs(timestep_adjustment_factor) != 1.0
        output_steps = model.timestep.parameters.output_steps
        nskip = output_steps.nskip.value
        output_steps.nskip.value = floor(Int, nskip/timestep_adjustment_factor)
    end
    return nothing
end

function get_timestep_adjustment_factor(model::ModelData)::Float64
    ivelocity_step_counter = model.bcs.parameters.velocity_step.ivelocity_step_counter.value
    velocity_step = model.bcs.parameters.velocity_step

    if ivelocity_step_counter == 0
        return velocity_step.timestep_adjustment_factor.value
    elseif ivelocity_step_counter == 1
        return velocity_step.timestep_second_adjustment_factor.value
    else
        return 1.0
    end
end

function update_ivelocity_step_counter!(model::ModelData)::Nothing
    model.bcs.parameters.velocity_step.ivelocity_step_counter.value += 1
    return nothing
end

end # module 