module VelocityFromStrainRate

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: load_parameters!
import EarthBox.InitializationTools: update_use_flags

const PDATA = get_eb_parameters()

struct ValidInputNames
    iuse_strain_rate::Symbol
    strain_rate_bc::Symbol
end

""" 
    initialize!(model::ModelData; kwargs...)::Nothing

Initialize velocity from strain rate boundary condition parameters.

# Arguments
- `model::`[`ModelData`](@ref ModelData): Model data object containing model parameters and arrays.

# Optional Keyword Arguments
- `iuse_strain_rate::Union{Bool, Nothing}`: 
    - $(PDATA.iuse_strain_rate.description)
- `strain_rate_bc::Union{Float64, Nothing}`: 
    - $(PDATA.strain_rate_bc.description)

# Returns
- `Nothing`

"""
function initialize!(
    model::ModelData;
    kwargs...
)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    if model.bcs.parameters.velocity.iuse_strain_rate.value == 1
        full_velocity_extension = calculate_velocity_from_strain_rate!(model)
        set_full_extension_velocity!(model, full_velocity_extension)
    end
    return nothing
end

function calculate_velocity_from_strain_rate!(model::ModelData)::Float64
    strain_rate = model.bcs.parameters.velocity.strain_rate_bc.value
    model_width = model.grids.parameters.geometry.xsize.value
    return strain_rate * model_width
end

function set_full_extension_velocity!(
    model::ModelData,
    full_velocity_extension::Float64
)::Nothing
    model.bcs.parameters.velocity.full_velocity_extension.value = full_velocity_extension
    return nothing
end

end # module 