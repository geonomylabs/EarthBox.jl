module Uniform

import EarthBox.ModelDataContainer: ModelData

"""
    initialize!(model::ModelData)

Calculate initial uniform temperature for each marker.

# Arguments
- `model::ModelData`: The model data structure containing marker information

# Updated Arrays
- `model.markers.arrays.thermal.marker_TK.array`: Temperature of markers in Kelvins
"""
function initialize!(model::ModelData)
    temperature_uniform = model.heat_equation.parameters.initial_condition.temperature_uniform.value
    marker_temperature = model.markers.arrays.thermal.marker_TK.array
    mxnum = model.markers.parameters.distribution.mxnum.value
    mynum = model.markers.parameters.distribution.mynum.value
    
    imarker = 1  # Julia uses 1-based indexing
    for _ixm in 1:mxnum
        for _iym in 1:mynum
            marker_temperature[imarker] = temperature_uniform
            imarker += 1
        end
    end
end

end # module 