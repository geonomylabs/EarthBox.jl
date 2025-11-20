module Plume

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: load_parameters!

const PDATA = get_eb_parameters()

struct ValidInputNames
    plume_radius::Symbol
    plume_center_x::Symbol
    plume_center_y::Symbol
    plume_head_thick::Symbol
    delta_temperature_plume::Symbol
    iuse_plume::Symbol
end

""" Initialize plume geometry.

# Keyword arguments:
- `plume_radius::Float64`:
    - $(PDATA.plume_radius.description)
- `plume_center_x::Float64`:
    - $(PDATA.plume_center_x.description)
- `plume_center_y::Float64`:
    - $(PDATA.plume_center_y.description)
- `plume_head_thick::Float64`:
    - $(PDATA.plume_head_thick.description)
- `delta_temperature_plume::Union{Float64, Nothing}`:
    - $(PDATA.delta_temperature_plume.description)
"""
function initialize!(model::ModelData; kwargs...)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    # Turn on plume when initialization geometry is used
    model.geometry.parameters.plume.iuse_plume.value = 1
    return nothing
end

end 