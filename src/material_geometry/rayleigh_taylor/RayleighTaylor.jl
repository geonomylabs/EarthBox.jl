module RayleighTaylor

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: load_parameters!

const PDATA = get_eb_parameters()

struct ValidInputNames
    depth_interface_h1::Symbol
    wave_length_lambda::Symbol
    amplitude_initial::Symbol
end

""" Initialize Rayleigh Taylor geometry.

# Keyword arguments:
- `depth_interface_h1::Float64`:
    - $(PDATA.depth_interface_h1.description)
- `wave_length_lambda::Float64`:
    - $(PDATA.wave_length_lambda.description)
- `amplitude_initial::Float64`:
    - $(PDATA.amplitude_initial.description)
"""
function initialize!(model::ModelData; kwargs...)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    return nothing
end

end 