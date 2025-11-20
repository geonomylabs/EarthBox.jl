module HotBox

import EarthBox.ModelDataContainer: ModelData
import EarthBox.ConversionFuncs: celsius_to_kelvin

function initialize!(model::ModelData)
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marker_temperature = model.markers.arrays.thermal.marker_TK.array

    temperature_top = model.bcs.parameters.temperature.temperature_top.value
    temperature_bottom = model.bcs.parameters.temperature.temperature_bottom.value

    y_sealevel = model.topography.parameters.sealevel.y_sealevel.value
    ysize = model.grids.parameters.geometry.ysize.value

    temperature_box = 1976.0
    xmin_box = 125_000.0
    ymin_box = 30_000.0
    width_box = 40_000.0
    mxnum = model.markers.parameters.distribution.mxnum.value
    mynum = model.markers.parameters.distribution.mynum.value
    imarker = 1  # Julia uses 1-based indexing
    for _ixm in 1:mxnum
        for _iym in 1:mynum
            x_marker = marker_x[imarker]
            y_marker = marker_y[imarker]
            temperature_kelvins = calculate_temperature(
                x_marker, y_marker, xmin_box, ymin_box, width_box,
                temperature_box, temperature_top, temperature_bottom,
                y_sealevel, ysize
            )
            marker_temperature[imarker] = temperature_kelvins
            imarker += 1
        end
    end
end

"""
    calculate_temperature(
        x_marker::Float64,
        y_marker::Float64,
        xmin_box::Float64,
        ymin_box::Float64,
        width_box::Float64,
        temperature_box::Float64,
        temperature_top::Float64,
        temperature_bottom::Float64,
        y_sealevel::Float64,
        ysize::Float64
    )::Float64

Calculate marker temperature for the hot box benchmark.
"""
function calculate_temperature(
    x_marker::Float64,
    y_marker::Float64,
    xmin_box::Float64,
    ymin_box::Float64,
    width_box::Float64,
    temperature_box::Float64,
    temperature_top::Float64,
    temperature_bottom::Float64,
    y_sealevel::Float64,
    ysize::Float64
)::Float64
    dT_dy = (temperature_bottom - temperature_top) / (ysize - y_sealevel)
    temperature = temperature_top + dT_dy * (y_marker - y_sealevel)
    if in_box(x_marker, y_marker, xmin_box, ymin_box, width_box)
        temperature = temperature_box
    end
    return temperature
end

function in_box(
    x_marker::Float64,
    y_marker::Float64,
    xmin_box::Float64,
    ymin_box::Float64,
    width_box::Float64
)::Bool
    xmax_box = xmin_box + width_box
    ymax_box = ymin_box + width_box
    return xmin_box <= x_marker <= xmax_box && ymin_box <= y_marker <= ymax_box
end

end # module 