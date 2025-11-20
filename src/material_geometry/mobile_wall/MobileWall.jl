module MobileWall

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: load_parameters!

const PDATA = get_eb_parameters()

struct ValidInputNames
    x_left_mobile_wall::Symbol
    x_right_mobile_wall::Symbol
    y_top_mobile_wall::Symbol
    y_bottom_mobile_wall::Symbol
    plate_extension_width::Symbol
    plate_extension_thickness::Symbol
end

""" Initialize mobile wall geometry.

# Keyword arguments:
- `x_left_mobile_wall::Union{Float64, Nothing}`:
    - $(PDATA.x_left_mobile_wall.description)
- `x_right_mobile_wall::Union{Float64, Nothing}`:
    - $(PDATA.x_right_mobile_wall.description)
- `y_top_mobile_wall::Union{Float64, Nothing}`:
    - $(PDATA.y_top_mobile_wall.description)
- `y_bottom_mobile_wall::Union{Float64, Nothing}`:
    - $(PDATA.y_bottom_mobile_wall.description)
- `plate_extension_width::Union{Float64, Nothing}`:
    - $(PDATA.plate_extension_width.description)
- `plate_extension_thickness::Union{Float64, Nothing}`:
    - $(PDATA.plate_extension_thickness.description)
"""
function initialize!(model::ModelData; kwargs...)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    return nothing
end

end 