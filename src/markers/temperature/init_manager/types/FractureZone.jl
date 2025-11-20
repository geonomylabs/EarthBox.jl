module FractureZone

using EarthBox.ModelDataContainer: ModelData
import EarthBox.ConversionFuncs: millions_of_years_to_seconds
using SpecialFunctions: erf

"""
    struct ThermalParameters

Thermal parameters for the fracture zone model.

# Fields
- `temperature_top::Float64`: Temperature at the top in Kelvin
- `temperature_bottom::Float64`: Temperature at the bottom in Kelvin
- `adiabatic_gradient::Float64`: Adiabatic gradient in Kelvin per kilometer
- `x_fracture_zone_end::Float64`: X-coordinate of fracture zone end in meters
- `age_left::Float64`: Age of left lithosphere in seconds
- `age_right::Float64`: Age of right lithosphere in seconds
- `sticky_thickness::Float64`: Thickness of sticky air layer in meters
- `thermal_lithosphere_depth_left::Float64`: Thermal lithosphere depth on left in meters
- `thermal_lithosphere_depth_right::Float64`: Thermal lithosphere depth on right in meters
- `thermal_diffusivity::Float64`: Thermal diffusivity in square meters per second
- `ysize::Float64`: Model size in y-direction in meters
"""
struct ThermalParameters
    temperature_top::Float64  # K
    temperature_bottom::Float64  # K
    adiabatic_gradient::Float64  # K/km
    x_fracture_zone_end::Float64  # m
    age_left::Float64  # s
    age_right::Float64  # s
    sticky_thickness::Float64  # m
    thermal_lithosphere_depth_left::Float64  # m
    thermal_lithosphere_depth_right::Float64  # m
    thermal_diffusivity::Float64  # m^2/s
    ysize::Float64  # model size in y-direction m
end

function initialize!(model::ModelData)::Nothing
    thermal_parameters = get_thermal_parameters(model)

    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marker_temperature = model.markers.arrays.thermal.marker_TK.array
    mxnum = model.markers.parameters.distribution.mxnum.value
    mynum = model.markers.parameters.distribution.mynum.value

    imarker = 1  # Julia uses 1-based indexing
    for _ixm in 1:mxnum
        for _iym in 1:mynum
            x_marker = marker_x[imarker]
            y_marker = marker_y[imarker]
            marker_temperature[imarker] = calculate_temperature(
                y_marker, x_marker, thermal_parameters
            )
            imarker += 1
        end
    end
    return nothing
end

function get_thermal_parameters(model::ModelData)::ThermalParameters
    init_params = model.heat_equation.parameters.initial_condition

    age_left = init_params.age_lithosphere_left.value
    age_left = millions_of_years_to_seconds(age_left)

    age_right = init_params.age_lithosphere_right.value
    age_right = millions_of_years_to_seconds(age_right)

    temp_params = model.bcs.parameters.temperature
    geom_params = model.geometry.parameters
    thermal_parameters = ThermalParameters(
        temp_params.temperature_top.value,
        temp_params.temperature_bottom.value,
        init_params.adiabatic_gradient.value,
        geom_params.fracture_zone.x_fracture_zone_end.value,
        age_left,
        age_right,
        geom_params.sticky_air_geometry.thick_air.value,
        init_params.thermal_lithosphere_depth_left.value,
        init_params.thermal_lithosphere_depth_right.value,
        init_params.thermal_diffusivity.value,
        model.grids.parameters.geometry.ysize.value
    )

    return thermal_parameters
end

"""
    calculate_temperature(
        y_marker::Float64,
        x_marker::Float64,
        thermal_parameters::ThermalParameters
    )::Float64

Calculate temperature using a simple fracture zone model.

# Arguments
- `y_marker::Float64`: Y-coordinate of marker in meters
- `x_marker::Float64`: X-coordinate of marker in meters
- `thermal_parameters::ThermalParameters`: Thermal parameters for the fracture zone model

# Returns
- `Float64`: Temperature in Kelvin at the marker location
"""
function calculate_temperature(
    y_marker::Float64,
    x_marker::Float64,
    thermal_parameters::ThermalParameters
)::Float64
    kappa = thermal_parameters.thermal_diffusivity
    temperature_top = thermal_parameters.temperature_top
    temperature_bottom = thermal_parameters.temperature_bottom
    adiabatic_gradient = thermal_parameters.adiabatic_gradient
    x_fracture_zone_end = thermal_parameters.x_fracture_zone_end
    age_left = thermal_parameters.age_left
    age_right = thermal_parameters.age_right
    sticky_thickness = thermal_parameters.sticky_thickness
    thermal_lithosphere_depth_left = thermal_parameters.thermal_lithosphere_depth_left
    thermal_lithosphere_depth_right = thermal_parameters.thermal_lithosphere_depth_right
    ysize = thermal_parameters.ysize

    # Depth below sticky air/water layer
    sub_sticky_depth = y_marker - sticky_thickness
    # Asthenosphere thickness
    left_asthenosphere_thickness = ysize - thermal_lithosphere_depth_left
    right_asthenosphere_thickness = ysize - thermal_lithosphere_depth_right
    # Initialize with adiabatic temperature gradient
    temperature = temperature_bottom - (
        adiabatic_gradient*(ysize - y_marker)/1000.0)
    # Sticky water condition
    if y_marker <= sticky_thickness
        temperature = temperature_top
    # Oceanic geotherm to the left of the fracture zone
    elseif (
        0 < x_marker < x_fracture_zone_end
        && sticky_thickness < y_marker < thermal_lithosphere_depth_left
    )
        temperature = calc_lithosphere_temperature(
            temperature_top, temperature_bottom, adiabatic_gradient,
            left_asthenosphere_thickness, kappa, sub_sticky_depth, age_left
        )
    # Oceanic geotherm to the right of the fracture zone
    elseif (
        x_marker > x_fracture_zone_end
        && sticky_thickness < y_marker < thermal_lithosphere_depth_right
    )
        temperature = calc_lithosphere_temperature(
            temperature_top, temperature_bottom, adiabatic_gradient,
            right_asthenosphere_thickness, kappa, sub_sticky_depth, age_right
        )
    end

    return temperature
end

"""
    calc_lithosphere_temperature(
        temperature_top::Float64,
        temperature_bottom::Float64,
        adiabatic_gradient::Float64,
        asthenosphere_thickness::Float64,
        kappa::Float64,
        depth::Float64,
        age::Float64
    )::Float64

Calculate temperature for a simple lithosphere geotherm.

# Arguments
- `temperature_top::Float64`: Surface temperature in Kelvin
- `temperature_bottom::Float64`: Bottom temperature in Kelvin
- `adiabatic_gradient::Float64`: Adiabatic gradient in Kelvin per kilometer
- `asthenosphere_thickness::Float64`: Thickness of the asthenosphere in meters
- `kappa::Float64`: Thermal diffusivity in square meters per second
- `depth::Float64`: Depth in meters
- `age::Float64`: Age in seconds

# Returns
- `Float64`: Temperature in Kelvin
"""
function calc_lithosphere_temperature(
    temperature_top::Float64,
    temperature_bottom::Float64,
    adiabatic_gradient::Float64,
    asthenosphere_thickness::Float64,
    kappa::Float64,
    depth::Float64,
    age::Float64
)::Float64
    temperature_base_lith = (
        temperature_bottom
        - adiabatic_gradient*asthenosphere_thickness/1000.0
    )
    temperature = (
        temperature_base_lith + (temperature_top - temperature_base_lith)
        *(1 - erf((depth)/(2*(kappa*age)^0.5)))
    )
    return temperature
end

end # module 