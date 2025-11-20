module ReferenceLithosphere

import EarthBox.ConversionFuncs: celsius_to_kelvin
import EarthBox.Markers.MarkerTemperature.InitManager.TempIniStructs: 
    LayerThickness, Perturbation, TemperatureBCs, Conductivity, HeatProduction, InternalTemperature
import EarthBox.Markers.MarkerTemperature.InitManager.FourLinearSegments: 
    calculate_temperature as calculate_temperature_four_linear_segments
import EarthBox.Markers.MarkerTemperature.InitManager.AnalyticalThreeLayer: 
    calculate_temperature as calculate_temperature_analytical_three_layer
import EarthBox.RockProperties.DensityModel.UpdateManager.DensityLiao14: marker_density_model
import ..DensityProps: get_density_model_based_on_inputs
import ..DensityProps: get_density_props_at_y_depth
import ..DensityProps: DensityModel

struct LithosphereThicknesses
    total_column_thickness_meters::Float64  # Total thickness of the lithosphere column in meters
    thickness_upr_cont_crust_meters::Float64  # Thickness of the upper continental crust in meters
    thickness_lwr_cont_crust_meters::Float64  # Thickness of the lower continental crust in meters
    thickness_lithosphere_meters::Float64  # Thickness of the lithosphere in meters
    dy_meters::Float64  # Grid spacing in meters
end

struct LithosphereThermalProps
    temperature_top_celsius::Float64  # Temperature at the top of the reference column in Celsius
    temperature_moho_celsius::Float64  # Temperature at the Moho in Celsius
    temperature_base_lith_celsius::Float64  # Temperature at the base of the lithosphere in Celsius
    adiabatic_gradient_kelvin_km::Float64  # Adiabatic gradient in Kelvin per kilometer
    conductivity_upper_crust::Float64  # Thermal conductivity of the upper crust in W/m/K
    conductivity_lower_crust::Float64  # Thermal conductivity of the lower crust in W/m/K
    conductivity_mantle::Float64  # Thermal conductivity of the mantle in W/m/K
    heat_production_upper_crust::Float64  # Heat production in the upper crust in W/m^3
    heat_production_lower_crust::Float64  # Heat production in the lower crust in W/m^3
    heat_production_mantle::Float64  # Heat production in the mantle in W/m^3
    thickness_thermal_lithosphere::Float64  # Thickness of the thermal lithosphere in meters
end

""" Make lithosphere model from user inputs using reference lithosphere.

# Keyword Arguments
- `iuse_linear_segments::Bool`: Use a temperature profile with four linear segments. If false, 
  an analytical three-layer temperature profile is used instead.
- `thickness_upr_cont_crust_meters::Float64`: Thickness of upper continental crust in meters.
- `thickness_lwr_cont_crust_meters::Float64`: Thickness of lower continental crust in meters.
- `thickness_lithosphere_meters::Float64`: Thickness of lithosphere in meters.
- `thickness_thermal_lithosphere::Float64`: Thickness of thermal lithosphere in meters.
- `thickness_asthenosphere_meters::Float64`: Thickness of asthenosphere in meters.
- `dy_meters::Float64`: Grid spacing in meters.
- `expansivity::Float64`: Thermal expansivity in 1/K.
- `compressibility::Float64`: Compressibility in 1/Pa.
- `density_upper_continental_crust::Float64`: Density of upper continental crust in kg/m^3.
- `density_lower_continental_crust::Float64`: Density of lower continental crust in kg/m^3.
- `density_mantle_lithosphere::Float64`: Density of mantle lithosphere in kg/m^3.
- `density_asthenosphere::Float64`: Density of asthenosphere in kg/m^3.
- `temperature_top_celsius::Float64`: Temperature at top of lithosphere in Celsius.
- `temperature_moho_celsius::Float64`: Temperature at Moho in Celsius.
- `temperature_base_lith_celsius::Float64`: Temperature at base of lithosphere in Celsius.
- `adiabatic_gradient_kelvin_km::Float64`: Adiabatic gradient in Kelvin/km.
- `conductivity_upper_crust::Float64`: Thermal conductivity of upper crust in W/m/K.
- `conductivity_lower_crust::Float64`: Thermal conductivity of lower crust in W/m/K.
- `conductivity_mantle::Float64`: Thermal conductivity of mantle in W/m/K.
- `heat_production_upper_crust::Float64`: Heat production in upper crust in W/m^3.
- `heat_production_lower_crust::Float64`: Heat production in lower crust in W/m^3.
- `heat_production_mantle::Float64`: Heat production in mantle in W/m^3.

# Returns
- `gridy::Vector{Float64}`: Grid in y-direction.
- `temp_gridy::Vector{Float64}`: Temperature (Kelvin) grid in y-direction.
- `density_gridy::Vector{Float64}`: Density (kg/m^3) grid in y-direction.
- `pressure_gridy::Vector{Float64}`: Pressure (Pascals) grid in y-direction.
- `lith_thicknesses::LithosphereThicknesses`: Lithosphere thicknesses.
"""
function make_lithosphere_model_from_user_inputs(;
    kwargs...
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, LithosphereThicknesses}
    
    iuse_linear_segments = get(kwargs, :iuse_linear_segments, true)
    thickness_upr_cont_crust_meters = get(kwargs, :thickness_upr_cont_crust_meters, 22_000.0)
    thickness_lwr_cont_crust_meters = get(kwargs, :thickness_lwr_cont_crust_meters, 10_000.0)
    thickness_lithosphere_meters = get(kwargs, :thickness_lithosphere_meters, 125_000.0)
    thickness_thermal_lithosphere = get(kwargs, :thickness_thermal_lithosphere, 120_000.0)
    thickness_asthenosphere_meters = get(kwargs, :thickness_asthenosphere_meters, 20_000.0)
    dy_meters = get(kwargs, :dy_meters, 100.0)
    expansivity = get(kwargs, :expansivity, 3.0e-5)
    compressibility = get(kwargs, :compressibility, 1.0e-11)
    density_upper_continental_crust = get(kwargs, :density_upper_continental_crust, 2800.0)
    density_lower_continental_crust = get(kwargs, :density_lower_continental_crust, 2800.0)
    density_mantle_lithosphere = get(kwargs, :density_mantle_lithosphere, 2800.0)
    density_asthenosphere = get(kwargs, :density_asthenosphere, 3300.0)
    temperature_top_celsius = get(kwargs, :temperature_top_celsius, 0.0)
    temperature_moho_celsius = get(kwargs, :temperature_moho_celsius, 600.0)
    temperature_base_lith_celsius = get(kwargs, :temperature_base_lith_celsius, 1330.0)
    adiabatic_gradient_kelvin_km = get(kwargs, :adiabatic_gradient_kelvin_km, 0.4)
    conductivity_upper_crust = get(kwargs, :conductivity_upper_crust, 2.5)
    conductivity_lower_crust = get(kwargs, :conductivity_lower_crust, 2.5)
    conductivity_mantle = get(kwargs, :conductivity_mantle, 2.25)
    heat_production_upper_crust = get(kwargs, :heat_production_upper_crust, 9e-7)
    heat_production_lower_crust = get(kwargs, :heat_production_lower_crust, 9e-7)
    heat_production_mantle = get(kwargs, :heat_production_mantle, 0.0)

    density_model = get_density_model_based_on_inputs(
        expansivity=expansivity,
        compressibility=compressibility,
        density_upper_continental_crust=density_upper_continental_crust,
        density_lower_continental_crust=density_lower_continental_crust,
        density_mantle_lithosphere=density_mantle_lithosphere,
        density_asthenosphere=density_asthenosphere
    )

    total_column_thickness_meters = (
        thickness_upr_cont_crust_meters +
        thickness_lwr_cont_crust_meters +
        thickness_lithosphere_meters +
        thickness_asthenosphere_meters
    )

    lith_thicknesses = LithosphereThicknesses(
        total_column_thickness_meters,
        thickness_upr_cont_crust_meters,
        thickness_lwr_cont_crust_meters,
        thickness_lithosphere_meters,
        dy_meters
    )

    lith_thermal = LithosphereThermalProps(
        temperature_top_celsius,
        temperature_moho_celsius,
        temperature_base_lith_celsius,
        adiabatic_gradient_kelvin_km,
        conductivity_upper_crust,
        conductivity_lower_crust,
        conductivity_mantle,
        heat_production_upper_crust,
        heat_production_lower_crust,
        heat_production_mantle,
        thickness_thermal_lithosphere
    )

    gridy, temp_gridy, density_gridy, pressure_gridy = reference_lithosphere(
        density_model,
        lith_thicknesses,
        lith_thermal,
        number_of_pressure_iterations=3,
        iuse_linear_segments=iuse_linear_segments,
        gravity=9.8
    )

    return gridy, temp_gridy, density_gridy, pressure_gridy, lith_thicknesses
end

""" Define the reference lithosphere.

# Arguments
- `density_model::DensityModel`: Density model.
- `lith_thicknesses::LithosphereThicknesses`: Reference lithosphere thicknesses.
- `lith_thermal::LithosphereThermalProps`: Reference lithosphere thermal properties.
- `number_of_pressure_iterations::Int`: Number of pressure iterations.
- `iuse_linear_segments::Int64`: Use a temperature profile with four linear segments.
- `gravity::Float64`: Gravitational acceleration in m/s^2.

# Returns
- `gridy::Vector{Float64}`: Array of y-coordinates (meters).
- `temp_gridy::Vector{Float64}`: Array of temperatures (Kelvin) at gridy nodes.
- `density_gridy::Vector{Float64}`: Array of densities (kg/m^3) at gridy nodes.
- `pressure_gridy::Vector{Float64}`: Array of pressures (Pascals) at gridy nodes.
"""
function reference_lithosphere(
    density_model::DensityModel,
    lith_thicknesses::LithosphereThicknesses,
    lith_thermal::LithosphereThermalProps;
    number_of_pressure_iterations::Int = 3,
    iuse_linear_segments::Int64 = 1,
    gravity::Float64 = 9.8
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}

    total_column_thickness_meters = lith_thicknesses.total_column_thickness_meters
    thickness_upr_cont_crust_meters = lith_thicknesses.thickness_upr_cont_crust_meters
    thickness_lwr_cont_crust_meters = lith_thicknesses.thickness_lwr_cont_crust_meters
    thickness_lithosphere_meters = lith_thicknesses.thickness_lithosphere_meters
    dy_meters = lith_thicknesses.dy_meters

    temperature_top_celsius = lith_thermal.temperature_top_celsius
    temperature_moho_celsius = lith_thermal.temperature_moho_celsius
    temperature_base_lith_celsius = lith_thermal.temperature_base_lith_celsius
    adiabatic_gradient_kelvin_km = lith_thermal.adiabatic_gradient_kelvin_km
    conductivity_upper_crust = lith_thermal.conductivity_upper_crust
    conductivity_lower_crust = lith_thermal.conductivity_lower_crust
    conductivity_mantle = lith_thermal.conductivity_mantle
    heat_production_upper_crust = lith_thermal.heat_production_upper_crust
    heat_production_lower_crust = lith_thermal.heat_production_lower_crust
    heat_production_mantle = lith_thermal.heat_production_mantle
    thickness_thermal_lithosphere = lith_thermal.thickness_thermal_lithosphere

    thickness_mantle_lithosphere = (
        thickness_lithosphere_meters -
        (thickness_upr_cont_crust_meters + thickness_lwr_cont_crust_meters)
    )

    gridy = make_reference_lithosphere_grid(total_column_thickness_meters, dy_meters)

    layer_thickness = LayerThickness(
        0.0,
        thickness_upr_cont_crust_meters,
        thickness_lwr_cont_crust_meters,
        thickness_mantle_lithosphere/3.0,
        thickness_mantle_lithosphere/3.0,
        thickness_mantle_lithosphere/3.0,
        thickness_thermal_lithosphere
    )

    conductivity = Conductivity(
        conductivity_upper_crust,
        conductivity_lower_crust,
        conductivity_mantle
    )

    heat_production = HeatProduction(
        heat_production_upper_crust,
        heat_production_lower_crust,
        heat_production_mantle
    )

    perturbation = Perturbation(0.0, 0.0)

    temperature_bottom_kelvin = calculate_temperature_bottom(
        total_column_thickness_meters,
        adiabatic_gradient_kelvin_km,
        thickness_lithosphere_meters,
        celsius_to_kelvin(temperature_base_lith_celsius)
    )

    temperature_bcs = TemperatureBCs(
        celsius_to_kelvin(temperature_top_celsius),
        temperature_bottom_kelvin,
        celsius_to_kelvin(temperature_base_lith_celsius)
    )

    internal_temperature = InternalTemperature(
        celsius_to_kelvin(temperature_top_celsius),
        celsius_to_kelvin(temperature_moho_celsius)
    )

    temp_gridy = make_temperature_grid(
        gridy,
        temperature_bcs,
        internal_temperature,
        layer_thickness,
        perturbation,
        conductivity,
        heat_production,
        iuse_linear_segments=iuse_linear_segments
    )

    # Initialize pressure grid assuming zero pressure
    pressure_gridy = Vector{Float64}(undef, length(gridy)) #zeros(Float64, length(gridy))
    density_gridy = Vector{Float64}(undef, length(gridy)) #zeros(Float64, length(gridy))
    Threads.@threads for i in 1:length(gridy)
        pressure_gridy[i] = 0.0
        density_gridy[i] = 0.0
    end

    for _ in 1:number_of_pressure_iterations
        density_gridy = calculate_density_grid(
            gridy,
            layer_thickness,
            density_model,
            temp_gridy,
            pressure_gridy
        )

        pressure_gridy = calculate_pressure_grid(gridy, density_gridy, gravity=gravity)
    end

    return gridy, temp_gridy, density_gridy, pressure_gridy
end

""" Make the reference lithosphere grid.

# Arguments
- `model_thickness_meters::Float64`: Thickness of the model in meters.
- `dy_meters::Float64`: Grid spacing in meters.

# Returns
- `gridy::Vector{Float64}`: Array of y-coordinates (meters).
"""
function make_reference_lithosphere_grid(
    model_thickness_meters::Float64,
    dy_meters::Float64
)::Vector{Float64}
    nnodes = floor(Int, model_thickness_meters/dy_meters) + 1
    gridy = Vector{Float64}(undef, nnodes)
    for i in 1:nnodes
        gridy[i] = (i-1) * dy_meters
    end
    return gridy
end

""" Make the temperature grid.

# Arguments
- `gridy::Vector{Float64}`: Array of y-coordinates (meters).
- `temperature_bcs::TemperatureBCs`: Boundary conditions for temperature.
- `internal_temperature::InternalTemperature`: Internal temperature struct.
- `layer_thickness::LayerThickness`: Layer thicknesses.
- `perturbation::Perturbation`: Perturbation parameters.
- `conductivity::Conductivity`: Conductivity parameters.
- `heat_production::HeatProduction`: Heat production parameters.
- `iuse_linear_segments::Int64`: Use a temperature profile with four linear segments.
- `adiabatic_gradient_kelvin_km::Float64`: Adiabatic gradient in Kelvin per kilometer.

# Returns
- `temp_gridy::Vector{Float64}`: Array of temperatures (Kelvin) at gridy nodes.
"""
function make_temperature_grid(
    gridy::Vector{Float64},
    temperature_bcs::TemperatureBCs,
    internal_temperature::InternalTemperature,
    layer_thickness::LayerThickness,
    perturbation::Perturbation,
    conductivity::Conductivity,
    heat_production::HeatProduction;
    iuse_linear_segments::Int64 = 1,
    adiabatic_gradient_kelvin_km::Float64 = 0.4
)::Vector{Float64}
    ysize = gridy[end]
    x_location = 0.0
    nnodes = length(gridy)
    temp_gridy = Vector{Float64}(undef, nnodes)
    
    for i in 1:nnodes
        y_location = gridy[i]
        if iuse_linear_segments == 1
            temp_gridy[i] = calculate_temperature_four_linear_segments(
                x_location,
                y_location,
                1e6,
                ysize,
                temperature_bcs,
                internal_temperature,
                layer_thickness,
                perturbation
            )
        else
            temp_gridy[i] = calculate_temperature_analytical_three_layer(
                x_location,
                y_location,
                1e6,
                layer_thickness,
                conductivity,
                heat_production,
                temperature_bcs,
                perturbation,
                adiabatic_gradient_kelvin_km
            )
        end
    end

    return temp_gridy
end

""" Calculate the temperature at the bottom of the reference column.

# Arguments
- `reference_column_thickness_meters::Float64`: Thickness of the reference column in meters.
- `adiabatic_gradient_kelvin_km::Float64`: Adiabatic gradient in Kelvin per kilometer.
- `thickness_lithosphere_meters::Float64`: Thickness of the lithosphere in meters.
- `temperature_base_lith_kelvin::Float64`: Temperature at the base of the lithosphere in Kelvin.

# Returns
- `temperature_bottom_kelvin::Float64`: Temperature at the bottom of the reference column in Kelvin.
"""
function calculate_temperature_bottom(
    reference_column_thickness_meters::Float64,
    adiabatic_gradient_kelvin_km::Float64,
    thickness_lithosphere_meters::Float64,
    temperature_base_lith_kelvin::Float64
)::Float64
    thickness_asthenosphere_meters = (
        reference_column_thickness_meters - thickness_lithosphere_meters
    )

    adiabatic_gradient_kelvin_meters = adiabatic_gradient_kelvin_km/1000.0

    temperature_bottom_kelvin = (
        temperature_base_lith_kelvin +
        adiabatic_gradient_kelvin_meters * thickness_asthenosphere_meters
    )

    return temperature_bottom_kelvin
end

""" Calculate the pressure grid.

# Arguments
- `gridy::Vector{Float64}`: Array of y-coordinates (meters).
- `density_gridy::Vector{Float64}`: Array of densities (kg/m^3) at gridy nodes.
- `gravity::Float64`: Gravitational acceleration in m/s^2.

# Returns
- `pressure_gridy::Vector{Float64}`: Array of pressures (Pascals) at gridy nodes.
"""
function calculate_pressure_grid(
    gridy::Vector{Float64},
    density_gridy::Vector{Float64};
    gravity::Float64 = 9.8
)::Vector{Float64}
    nnodes = length(gridy)
    pressure_gridy = Vector{Float64}(undef, nnodes)
    pressure_gridy[1] = 0.0
    for i in 2:nnodes
        dy = gridy[i] - gridy[i-1]
        density_cell_above = (density_gridy[i-1] + density_gridy[i])/2.0
        pressure_gridy[i] = pressure_gridy[i-1] + density_cell_above * gravity * dy
    end
    return pressure_gridy
end

""" Calculate the density grid.

# Arguments
- `gridy::Vector{Float64}`: Array of y-coordinates (meters).
- `layer_thickness::LayerThickness`: Layer thicknesses.
- `density_model::DensityModel`: Density model.
- `temp_gridy::Vector{Float64}`: Array of temperatures (Kelvin) at gridy nodes.
- `pressure_gridy::Vector{Float64}`: Array of pressures (Pascals) at gridy nodes.

# Returns
- `density_gridy::Vector{Float64}`: Array of densities (kg/m^3) at gridy nodes.
"""
function calculate_density_grid(
    gridy::Vector{Float64},
    layer_thickness::LayerThickness,
    density_model::DensityModel,
    temp_gridy::Vector{Float64},
    pressure_gridy::Vector{Float64}
)::Vector{Float64}
    nnodes = length(gridy)
    density_gridy = Vector{Float64}(undef, nnodes) #zeros(Float64, nnodes)
    
    for i in 1:nnodes
        y_location = gridy[i]
        temperature_kelvin = temp_gridy[i]
        pressure_pascal = pressure_gridy[i]

        (
            standard_density, expansivity, compressibility
        ) = get_density_props_at_y_depth(layer_thickness, density_model, y_location)

        density_gridy[i] = marker_density_model(
            standard_density,
            expansivity,
            compressibility,
            pressure_pascal,
            temperature_kelvin
        )
    end

    return density_gridy
end

""" Get the base depths of the layers.

# Arguments
- `layer_thicknesses::LithosphereThicknesses`: Lithosphere thicknesses.

# Returns
- `y_base_upr_crust::Float64`: Base depth of upper crust in meters.
- `y_base_lwr_crust::Float64`: Base depth of lower crust in meters.
- `y_base_upr_mantle_lithosphere::Float64`: Base depth of upper mantle lithosphere in meters.
- `y_base_mid_mantle_lithosphere::Float64`: Base depth of middle mantle lithosphere in meters.
- `y_base_lwr_mantle_lithosphere::Float64`: Base depth of lower mantle lithosphere in meters.
"""
function get_lithosphere_layer_base_depths(
    layer_thicknesses::LithosphereThicknesses
)::Tuple{Float64, Float64, Float64, Float64, Float64}
    y_base_upr_crust = layer_thicknesses.thickness_upr_cont_crust_meters

    y_base_lwr_crust = (
        layer_thicknesses.thickness_upr_cont_crust_meters +
        layer_thicknesses.thickness_lwr_cont_crust_meters
    )

    thickness_mantle_lithosphere = (
        layer_thicknesses.thickness_lithosphere_meters -
        (layer_thicknesses.thickness_upr_cont_crust_meters +
         layer_thicknesses.thickness_lwr_cont_crust_meters)
    )

    y_base_upr_mantle_lithosphere = (
        y_base_lwr_crust + thickness_mantle_lithosphere/3.0
    )

    y_base_mid_mantle_lithosphere = (
        y_base_upr_mantle_lithosphere + thickness_mantle_lithosphere/3.0
    )

    y_base_lwr_mantle_lithosphere = (
        y_base_mid_mantle_lithosphere + thickness_mantle_lithosphere/3.0
    )

    return (
        y_base_upr_crust,
        y_base_lwr_crust,
        y_base_upr_mantle_lithosphere,
        y_base_mid_mantle_lithosphere,
        y_base_lwr_mantle_lithosphere
    )
end

end # module 