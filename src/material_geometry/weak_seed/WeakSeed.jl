module WeakSeed

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: load_parameters!

const PDATA = get_eb_parameters()

struct ValidInputNames
    w_seed::Symbol
    x_seed::Symbol
    y_seed::Symbol
end

""" Initialize weak seed geometry.

# Keyword arguments:
- `w_seed::Float64`:
    - $(PDATA.w_seed.description)
- `x_seed::Float64`:
    - $(PDATA.x_seed.description)
- `y_seed::Float64`:
    - $(PDATA.y_seed.description)
"""
function initialize!(model::ModelData; kwargs...)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    # Turn on weak seed when initialization geometry is used
    model.geometry.parameters.weak_seed.iuse_weak_seed.value = 1
    return nothing
end

end 