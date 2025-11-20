module MarkerStressLimits

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: load_parameters!

const PDATA = get_eb_parameters()

struct ValidInputNames
    powerlaw_stress_min::Symbol
    yield_stress_min::Symbol
    yield_stress_max::Symbol
end

"""
    initialize!(model::ModelData; kwargs...)::Nothing

Initialize marker stress limits for dislocation creep and plasticity models.

# Arguments
- `model::ModelData`: The model data container containing the model parameters and arrays.

# Keyword Arguments
- `powerlaw_stress_min::Union{Float64, Nothing}=nothing`:
    - $(PDATA.powerlaw_stress_min.description)
- `yield_stress_min::Union{Float64, Nothing}=nothing`:
    - $(PDATA.yield_stress_min.description)
- `yield_stress_max::Union{Float64, Nothing}=nothing`:
    - $(PDATA.yield_stress_max.description)

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

end # module MarkerStressLimits 