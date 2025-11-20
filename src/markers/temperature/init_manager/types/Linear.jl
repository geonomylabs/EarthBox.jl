module Linear

import EarthBox.ModelDataContainer: ModelData

"""
    initialize(model::ModelData)

Calculate initial temperature for each marker based using linear model.

# Updated Arrays
- `model.markers.arrays.thermal.marker_TK.array`: Temperature of markers in Kelvins

# Background
The initial temperature of each marker is calculated using a linear increase with depth. The temperature of each marker is calculated using the equation:

T = t_top + y*(t_bottom - t_top)/ysize

where:
    T = temperature of marker in Kelvins
    t_top = temperature at the top of the domain in Kelvins
    t_bottom = temperature at the bottom of the domain in Kelvins
    y = y-coordinate of marker
    ysize = y-size of the domain
"""
function initialize!(model::ModelData)
    marker_y = model.markers.arrays.location.marker_y.array
    marker_temperature = model.markers.arrays.thermal.marker_TK.array

    t_top = model.bcs.parameters.temperature.temperature_top.value
    t_bottom = model.bcs.parameters.temperature.temperature_bottom.value
    ysize = model.grids.parameters.geometry.ysize.value

    mxnum = model.markers.parameters.distribution.mxnum.value
    mynum = model.markers.parameters.distribution.mynum.value
    imarker = 1  # Julia uses 1-based indexing
    for _ixm in 1:mxnum
        for _iym in 1:mynum
            y_marker = marker_y[imarker]
            temperature_kelvins = calculate_temperature(
                y_marker, ysize, t_top, t_bottom
            )
            marker_temperature[imarker] = temperature_kelvins
            imarker += 1
        end
    end
end

"""
    calculate_temperature(
        y_marker::Float64,
        ysize::Float64,
        temperature_top::Float64,
        temperature_bottom::Float64
    )::Float64

Calculate marker temperature using a linear increase with depth.

This function is used in sandbox, channel flow and Couette flow experiments.
"""
function calculate_temperature(
    y_marker::Float64,
    ysize::Float64,
    temperature_top::Float64,
    temperature_bottom::Float64
)::Float64
    dtdy = (temperature_bottom - temperature_top) / ysize
    temperature_kelvins = temperature_top + y_marker * dtdy
    return temperature_kelvins
end

end # module 