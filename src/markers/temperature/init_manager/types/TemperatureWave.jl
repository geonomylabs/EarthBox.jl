module TemperatureWave

import EarthBox.ModelDataContainer: ModelData

"""
    initialize!(model::ModelData)

Calculate initial temperature for each marker.

# Arguments
- `model::ModelData`: The model data structure containing marker information

# Updated Arrays
- `model.markers.arrays.thermal.marker_TK.array`: Temperature of markers in Kelvins
"""
function initialize!(model::ModelData)
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marker_temperature = model.markers.arrays.thermal.marker_TK.array
    temperature_top = model.bcs.parameters.temperature.temperature_top.value
    temperature_bottom = model.bcs.parameters.temperature.temperature_bottom.value
    temperature_of_wave = model.heat_equation.parameters.initial_condition.temperature_of_wave.value
    xsize = model.grids.parameters.geometry.xsize.value
    ysize = model.grids.parameters.geometry.ysize.value
    mxnum = model.markers.parameters.distribution.mxnum.value
    mynum = model.markers.parameters.distribution.mynum.value
    
    imarker = 1  # Julia uses 1-based indexing
    for _ixm in 1:mxnum
        for _iym in 1:mynum
            x_marker = marker_x[imarker]
            y_marker = marker_y[imarker]
            temperature_kelvins = calculate_temperature(
                x_marker, y_marker, xsize, ysize, temperature_top,
                temperature_bottom, temperature_of_wave
            )
            marker_temperature[imarker] = temperature_kelvins
            imarker += 1
        end
    end
end

"""
    calculate_temperature(
        x_marker::Float64, y_marker::Float64, xsize::Float64, ysize::Float64,
        temperature_top::Float64, temperature_bottom::Float64,
        temperature_of_wave::Float64
    )::Float64

Calculate marker temperature for the solid body rotation benchmark.
"""
function calculate_temperature(
    x_marker::Float64,
    y_marker::Float64,
    xsize::Float64,
    ysize::Float64,
    temperature_top::Float64,
    temperature_bottom::Float64,
    temperature_of_wave::Float64
)::Float64
    return calc_initial_temperature_for_solid_body_rotation_benchmark(
        x_marker,
        y_marker,
        xsize,
        ysize,
        temperature_top,
        temperature_bottom,
        temperature_of_wave
    )
end

"""
    calc_initial_temperature_for_solid_body_rotation_benchmark(
        xcoor::Float64, ycoor::Float64,
        xsize::Float64, ysize::Float64,
        temperature_top::Float64,
        temperature_bottom::Float64,
        temperature_wave::Float64
    )::Float64

Calculate initial temperature of solid body rotation benchmark.

This initial condition is a square zone of elevated temperature. The
elevated zone then moves as a rotational wave.
"""
function calc_initial_temperature_for_solid_body_rotation_benchmark(
    xcoor::Float64,
    ycoor::Float64,
    xsize::Float64,
    ysize::Float64,
    temperature_top::Float64,
    temperature_bottom::Float64,
    temperature_wave::Float64
)::Float64
    vertical_temperature_gradient = (temperature_bottom - temperature_top) / ysize
    temperature_kelvins = temperature_top + ycoor * vertical_temperature_gradient

    dx_fac = xcoor / xsize
    dy_fac = ycoor / ysize
    if 0.4 ≤ dx_fac ≤ 0.6 && 0.1 ≤ dy_fac ≤ 0.3
        temperature_kelvins = temperature_wave
    end

    return temperature_kelvins
end

end # module 