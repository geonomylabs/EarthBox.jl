module CrustalHole

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: load_parameters!

const PDATA = get_eb_parameters()

struct ValidInputNames
    xhole_start::Symbol
    xhole_middle::Symbol
    xhole_end::Symbol
    xhole_depth::Symbol
end

""" Initialize crustal hole geometry.

# Keyword arguments:
- `xhole_start::Union{Float64, Nothing}`:
    - $(PDATA.xhole_start.description)
- `xhole_middle::Union{Float64, Nothing}`:
    - $(PDATA.xhole_middle.description)
- `xhole_end::Union{Float64, Nothing}`:
    - $(PDATA.xhole_end.description)
- `xhole_depth::Union{Float64, Nothing}`:
    - $(PDATA.xhole_depth.description)
"""
function initialize!(model::ModelData; kwargs...)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    return nothing
end

end