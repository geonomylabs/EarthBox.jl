module Sandbox

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: load_parameters!

const PDATA = get_eb_parameters()

struct ValidInputNames
    nsand_layers::Symbol
    y_sand_air_interface::Symbol
    y_top_microbeads::Symbol
    y_bottom_microbeads::Symbol
    x_left_ramp::Symbol
    x_right_ramp::Symbol
    ramp_dip_deg::Symbol
    pdms_layer_width::Symbol
    pdms_layer_thickness::Symbol
end

""" Initialize sandbox geometry.

# Keyword arguments:
- `nsand_layers::Int`:
    - $(PDATA.nsand_layers.description)
- `y_sand_air_interface::Float64`:
    - $(PDATA.y_sand_air_interface.description)
- `y_top_microbeads::Union{Float64, Nothing}`:
    - $(PDATA.y_top_microbeads.description)
- `y_bottom_microbeads::Union{Float64, Nothing}`:
    - $(PDATA.y_bottom_microbeads.description)
- `x_left_ramp::Union{Float64, Nothing}`:
    - $(PDATA.x_left_ramp.description)
- `x_right_ramp::Union{Float64, Nothing}`:
    - $(PDATA.x_right_ramp.description)
- `ramp_dip_deg::Union{Float64, Nothing}`:
    - $(PDATA.ramp_dip_deg.description)
- `pdms_layer_width::Union{Float64, Nothing}`:
    - $(PDATA.pdms_layer_width.description)
- `pdms_layer_thickness::Union{Float64, Nothing}`:
    - $(PDATA.pdms_layer_thickness.description)
"""
function initialize!(model::ModelData; kwargs...)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    return nothing
end

end 