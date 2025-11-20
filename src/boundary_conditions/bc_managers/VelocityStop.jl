module VelocityStop

import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: load_parameters!
import EarthBox.InitializationTools: update_use_flags
import EarthBox.ConversionFuncs: seconds_to_years

struct ValidInputNames
    iuse_velocity_stop::Symbol
    velocity_stop_time::Symbol
end

""" 
    initialize!(model::ModelData; kwargs...)::Nothing

Initialize velocity stopping model.

# Arguments
- `model::`[`ModelData`](@ref ModelData): Model data object containing model parameters and arrays.

# Option Keyword Arguments
- `iuse_velocity_stop::Union{Bool, Nothing}`: 
    - Integer flag indicating whether to stop velocity boundary conditions at a defined time
- `velocity_stop_time::Union{Float64, Nothing}`: 
    - Time in million years when velocity boundary conditions are set to zero
"""
function initialize!(model::ModelData; kwargs...)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    return nothing
end

function reset_bc_for_veloc_stop!(model::ModelData)::Nothing
    iuse = model.bcs.parameters.velocity.iuse_velocity_stop.value
    if iuse == 1 && at_stop_time(model)
        println(">> Setting extensional/contractional velocity to zero.")
        update_extension_velocity!(model)
        update_contraction_velocity!(model)
        reset_ivelocity_stop_counter!(model)
    end
    return nothing
end

function at_stop_time(model::ModelData)::Bool
    timesum = model.timestep.parameters.main_time_loop.timesum.value
    velocity_stop_time = model.bcs.parameters.velocity.velocity_stop_time.value
    ivelocity_stop_counter = model.bcs.parameters.velocity.ivelocity_stop_counter.value
    timesum_myr = seconds_to_years(timesum)/1e6
    return timesum_myr >= velocity_stop_time && ivelocity_stop_counter == 0
end

function update_extension_velocity!(model::ModelData)::Nothing
    println(">> New full_velocity_extension was set to zero")
    model.bcs.parameters.velocity.full_velocity_extension.value = 0.0
    return nothing
end

function update_contraction_velocity!(model::ModelData)::Nothing
    println(">> New full_velocity_contraction was set to zero")
    model.bcs.parameters.velocity.full_velocity_contraction.value = 0.0
    return nothing
end

function reset_ivelocity_stop_counter!(model::ModelData)::Nothing
    # Integer counter for velocity jump. 1 means velocity stop has occurred.
    model.bcs.parameters.velocity.ivelocity_stop_counter.value = 1
    return nothing
end

end # module 