module HalfSpaceDoubleSpreading

using EarthBox.ModelDataContainer: ModelData
import EarthBox.ConversionFuncs: celsius_to_kelvin
using ..TempIniStructs
using ..GeothermStructs
using ..HalfSpaceCooling: calculate_temperature as calculate_half_space_temperature

function initialize!(model::ModelData)::Nothing
    thermal_diffusivity_m2_s = 1e-6

    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marker_temperature = model.markers.arrays.thermal.marker_TK.array

    temperature_top = model.bcs.parameters.temperature.temperature_top.value
    temperature_asthenosphere = celsius_to_kelvin(1330.0)

    xsize = model.grids.parameters.geometry.xsize.value
    x_ridge_left = xsize * 0.25
    x_ridge_right = xsize * 0.75
    xmid = xsize / 2.0

    sticky_thickness = model.geometry.parameters.sticky_air_geometry.thick_air.value

    half_spreading_rate = model.bcs.parameters.velocity.full_velocity_extension.value / 2.0

    mxnum = model.markers.parameters.distribution.mxnum.value
    mynum = model.markers.parameters.distribution.mynum.value

    imarker = 1  # Julia uses 1-based indexing
    for _ixm in 1:mxnum
        for _iym in 1:mynum
            x_marker = marker_x[imarker]
            y_marker = marker_y[imarker]
            
            # Determine which ridge to use based on x position
            x_ridge = x_marker <= xmid ? x_ridge_left : x_ridge_right
            
            temperature_kelvins = calculate_temperature(
                x_marker, y_marker, temperature_top, temperature_asthenosphere,
                x_ridge, sticky_thickness, half_spreading_rate,
                thermal_diffusivity_m2_s
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
        temperature_top::Float64,
        temperature_asthenosphere::Float64,
        x_ridge::Float64,
        sticky_thickness::Float64,
        half_spreading_rate::Float64,
        thermal_diffusivity_m2_s::Float64
    )::Float64

Calculate marker temperature for half-space spreading.

# Arguments
- `x_marker::Float64`: x-coordinate of marker in meters
- `y_marker::Float64`: y-coordinate of marker in meters
- `temperature_top::Float64`: Surface temperature in Kelvin
- `temperature_asthenosphere::Float64`: Asthenosphere temperature in Kelvin
- `x_ridge::Float64`: x-coordinate of ridge axis in meters
- `sticky_thickness::Float64`: Thickness of sticky air layer in meters
- `half_spreading_rate::Float64`: Half spreading rate in meters per second
- `thermal_diffusivity_m2_s::Float64`: Thermal diffusivity in square meters per second

# Returns
- `temperature_kelvins::Float64`: Temperature at the given depth below the seafloor in Kelvin
"""
function calculate_temperature(
    x_marker::Float64,
    y_marker::Float64,
    temperature_top::Float64,
    temperature_asthenosphere::Float64,
    x_ridge::Float64,
    sticky_thickness::Float64,
    half_spreading_rate::Float64,
    thermal_diffusivity_m2_s::Float64
)::Float64
    if y_marker <= sticky_thickness
        temperature_kelvins = temperature_top
    else
        depth_below_surface_meters = y_marker - sticky_thickness
        distance_from_ridge_axis_meters = abs(x_marker - x_ridge)
        temperature_kelvins = calculate_half_space_temperature(
            depth_below_surface_meters, temperature_top, temperature_asthenosphere,
            distance_from_ridge_axis_meters, half_spreading_rate,
            thermal_diffusivity_m2_s
        )
    end
    return temperature_kelvins
end

end # module 