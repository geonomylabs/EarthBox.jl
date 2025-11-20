module Velocity

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: load_parameters!
import EarthBox.ModelDataContainer: set_parameter!

const PDATA = get_eb_parameters()

struct ValidInputNames
    velocity::Symbol
    full_velocity_extension::Symbol
    full_velocity_extension_step1::Symbol
    full_velocity_extension_step2::Symbol
    full_velocity_contraction::Symbol
    velocity_shear::Symbol
    velocity_rotation::Symbol
    velocity_internal_x::Symbol
    velocity_internal_y::Symbol
    plate_thickness::Symbol
    smoothing_thickness::Symbol
end

""" 
    initialize!(model::ModelData; kwargs...)::Nothing

Initialize velocity boundary condition parameters.

# Arguments
- `model::`[`ModelData`](@ref ModelData): Model data object containing model parameters and arrays.

# Optional Keyword Arguments
- `velocity::Float64`: 
    - $(PDATA.velocity.description)
- `full_velocity_extension::Float64`: 
    - $(PDATA.full_velocity_extension.description)
- `full_velocity_extension_step1::Float64`: 
    - $(PDATA.full_velocity_extension_step1.description)
- `full_velocity_extension_step2::Float64`: 
    - $(PDATA.full_velocity_extension_step2.description)
- `full_velocity_contraction::Float64`: 
    - $(PDATA.full_velocity_contraction.description)
- `velocity_shear::Float64`: 
    - $(PDATA.velocity_shear.description)
- `velocity_rotation::Float64`: 
    - $(PDATA.velocity_rotation.description)
- `velocity_internal_x::Float64`: 
    - $(PDATA.velocity_internal_x.description)
- `velocity_internal_y::Float64`: 
    - $(PDATA.velocity_internal_y.description)
- `plate_thickness::Float64`: 
    - $(PDATA.plate_thickness.description)
- `smoothing_thickness::Float64`: 
    - $(PDATA.smoothing_thickness.description)

# Returns
- `Nothing`

"""
function initialize!(model::ModelData; kwargs...)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    update_velocity_parameters!(model, kwargs...)
    return nothing
end

function update_velocity_parameters!(model::ModelData, kwargs...)
    velocity = nothing
    for (key, val) in kwargs
        if key == :velocity
            velocity = val
            break
        end
    end
    if velocity !== nothing && typeof(velocity) == Float64
        model.bcs.parameters.velocity.full_velocity_extension.value = velocity
        model.bcs.parameters.velocity.full_velocity_contraction.value = velocity
        model.bcs.parameters.velocity.velocity_shear.value = velocity
        model.bcs.parameters.velocity.velocity_rotation.value = velocity
    end
    return nothing
end

end # module 