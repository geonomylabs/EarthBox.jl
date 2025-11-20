module StickyAirGeometry

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: load_parameters!

const PDATA = get_eb_parameters()

struct ValidInputNames
    thick_air::Symbol
end

""" Initialize sticky air geometry.

# Keyword arguments:
- `thick_air::Float64`:
    - $(PDATA.thick_air.description)
"""
function initialize!(model::ModelData; kwargs...)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    return nothing
end

end 