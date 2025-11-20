module BoxConvection

import EarthBox.ModelDataContainer: ModelData

"""
    initialize!(model::ModelData)::Nothing

Calculate initial temperature for each marker based on itype_temp.

# Arguments
- `model::ModelData`: Model data container containing marker and grid information

# Updates
- `model.markers.arrays.thermal.marker_TK.array`: Temperature of markers in Kelvins
"""
function initialize!(model::ModelData)::Nothing
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marker_temperature = model.markers.arrays.thermal.marker_TK.array
    
    xsize = model.grids.parameters.geometry.xsize.value
    ysize = model.grids.parameters.geometry.ysize.value
    
    temperature_bottom = model.bcs.parameters.temperature.temperature_bottom.value
    temperature_top = model.bcs.parameters.temperature.temperature_top.value
    
    mxnum = model.markers.parameters.distribution.mxnum.value
    mynum = model.markers.parameters.distribution.mynum.value
    
    imarker = 1  # Julia uses 1-based indexing
    for _ixm in 1:mxnum
        for _iym in 1:mynum
            x_marker = marker_x[imarker]
            y_marker = marker_y[imarker]
            temperature_kelvins = calculate_temperature(
                x_marker, y_marker, xsize, ysize, temperature_top,
                temperature_bottom
            )
            marker_temperature[imarker] = temperature_kelvins
            imarker += 1
        end
    end
    return nothing
end

"""
    calculate_temperature(
        x_marker::Float64,
        y_marker::Float64,
        xsize::Float64,
        ysize::Float64,
        temperature_top::Float64,
        temperature_bottom::Float64
    )::Float64

Calculate marker temperature for convection in a box benchmark.

Initial temperature profile with lateral perturbation.

# Arguments
- `x_marker::Float64`: X coordinate of marker
- `y_marker::Float64`: Y coordinate of marker
- `xsize::Float64`: Size of model in x direction
- `ysize::Float64`: Size of model in y direction
- `temperature_top::Float64`: Temperature at top boundary in Kelvin
- `temperature_bottom::Float64`: Temperature at bottom boundary in Kelvin

# Returns
- `Float64`: Temperature in Kelvin at the marker location
"""
function calculate_temperature(
    x_marker::Float64,
    y_marker::Float64,
    xsize::Float64,
    ysize::Float64,
    temperature_top::Float64,
    temperature_bottom::Float64
)::Float64
    xwt = x_marker/xsize
    ywt = y_marker/ysize

    ymid = 0.5
    if ywt > ymid
        ywt = (1.0 - ymid) + ymid*(abs(ywt - ymid)/(1.0 - ymid))^10
    else
        ywt = (1.0 - ymid) - (1.0 - ymid)*(abs(ywt - ymid)/ymid)^10
    end

    temperature_kelvins = (
        (
         temperature_top
         + (temperature_bottom - temperature_top)*ywt^(2.0 + xwt)
         + temperature_bottom
         + (temperature_top - temperature_bottom)*(1.0 - ywt)^(3.0 - xwt)
        )/2.0
    )
    return temperature_kelvins
end

end # module BoxConvection 