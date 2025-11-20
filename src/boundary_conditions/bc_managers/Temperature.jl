module Temperature

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: load_parameters!
import EarthBox.ModelDataContainer: set_parameter!

const PDATA = get_eb_parameters()

struct ValidInputNames
    temperature_top::Symbol
    temperature_bottom::Symbol
    temperature_left::Symbol
    temperature_right::Symbol
end

"""
    initialize!(model::ModelData; kwargs...)::Nothing

Initialize temperature boundary conditions.

# Arguments
- `model::`[`ModelData`](@ref ModelData): Model data object containing model parameters and arrays.

# Optional Keyword Arguments
- `temperature_top::Float64`: 
    - $(PDATA.temperature_top.description)
- `temperature_bottom::Float64`: 
    - $(PDATA.temperature_bottom.description)
- `temperature_left::Float64`: 
    - $(PDATA.temperature_left.description)
- `temperature_right::Float64`: 
    - $(PDATA.temperature_right.description)

# Returns
- `Nothing`

"""
function initialize!(
    model::ModelData;
    kwargs...
)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    return nothing
end

end # module 