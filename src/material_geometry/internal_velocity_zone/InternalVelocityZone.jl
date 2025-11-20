module InternalVelocityZone

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: load_parameters!

const PDATA = get_eb_parameters()

struct ValidInputNames
    x_vx_internal::Symbol
    y_min_vx_internal::Symbol
    y_max_vx_internal::Symbol
    x_vy_internal::Symbol
    y_min_vy_internal::Symbol
    y_max_vy_internal::Symbol
end

""" Initialize internal velocity zone geometry.

# Keyword arguments:
- x_vx_internal::Union{Float64, Nothing}:
    - $(PDATA.x_vx_internal.description)
- y_min_vx_internal::Union{Float64, Nothing}:
    - $(PDATA.y_min_vx_internal.description)
- y_max_vx_internal::Union{Float64, Nothing}:
    - $(PDATA.y_max_vx_internal.description)
- x_vy_internal::Union{Float64, Nothing}:
    - $(PDATA.x_vy_internal.description)
- y_min_vy_internal::Union{Float64, Nothing}:
    - $(PDATA.y_min_vy_internal.description)
- y_max_vy_internal::Union{Float64, Nothing}:
    - $(PDATA.y_max_vy_internal.description)
"""
function initialize!(model::ModelData; kwargs...)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    return nothing
end

end 