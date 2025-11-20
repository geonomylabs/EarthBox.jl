module Geotherm

import ..GeothermStructs: SteadyState3LayerData

"""
    linear_temp(z_depth, temperature_top, temperature_bottom, thick_thermal_lithosphere)

Linear geotherm.

Returns
-------
temperature::Float64
    Temperature.
"""
function linear_temp(
    z_depth::Float64,
    temperature_top::Float64,
    temperature_bottom::Float64,
    thick_thermal_lithosphere::Float64
)::Float64
    temperature = (
        temperature_top
        + (temperature_bottom - temperature_top)
            /thick_thermal_lithosphere*z_depth
    )
    return temperature
end

"""
    steady_state_3layer(z_depth::Float64, ss_3layer_data::SteadyState3LayerData)

Steady-state temperature for 3-layer model.

Analytical three layer geotherm model. Temperature boundary conditions are
defined at the top and bottom of the model. The temperature is calculated
at a given depth z_depth.

# Arguments
- `z_depth::Float64`
    - Sub-rock depth in meters.

- `ss_3layer_data::SteadyState3LayerData`
    - Struct for 3-layer model with the following attributes:
        - temperature_top: float : Temperature at the sticky-rock interface in
            Kelvin
        - temperature_bottom: float, Temperature at the bottom of the model in
            Kelvin
        - layer1_thickness: float, Thickness of the first layer in meters
        - layer2_thickness: float, Thickness of the second layer in meters
        - layer3_thickness: float, Thickness of the third layer in meters
        - thickness_total: float, Total thickness of the model in meters
        - layer1_heat_production: float, Heat production of the first layer in
            W/m^3
        - layer2_heat_production: float, Heat production of the second layer in
            W/m^3
        - layer3_heat_production: float, Heat production of the third layer in
            W/m^3
        - layer1_conductivity: float, Thermal conductivity of the first layer
            in W/mK
        - layer2_conductivity: float, Thermal conductivity of the second layer
            in W/mK
        - layer3_conductivity: float, Thermal conductivity of the third layer
            in W/mK
            named tuple for 3-layer model with parameters for temperature,
            thickness, heat production, and conductivity values for each layer.

# Returns
- `temperature::Float64`
    - Temperature in Kelvins
- `surface_heat_flow::Float64`
    - Surface heat flow in mW/m^2
"""
@inline function steady_state_3layer(
    z_depth::Float64,
    ss_3layer_data::SteadyState3LayerData
)::Tuple{Float64,Float64}
    temperature_top    = ss_3layer_data.temperature_top
    temperature_bottom = ss_3layer_data.temperature_bottom
    thick1             = ss_3layer_data.layer1_thickness
    thick2             = ss_3layer_data.layer2_thickness
    thick3             = ss_3layer_data.layer3_thickness
    thick_total        = ss_3layer_data.thickness_total
    hp1                = ss_3layer_data.layer1_heat_production
    hp2                = ss_3layer_data.layer2_heat_production
    hp3                = ss_3layer_data.layer3_heat_production
    k1                 = ss_3layer_data.layer1_conductivity
    k2                 = ss_3layer_data.layer2_conductivity
    k3                 = ss_3layer_data.layer3_conductivity

    @assert !isnan(thick1) "Layer 1 thickness (most likely upper crust) is NaN. Make sure you initialized the EarthLayering material geometry before initializing the temperature."
    @assert !isnan(thick2) "Layer 2 thickness (most likely lower crust) is NaN. Make sure you initialized the EarthLayering material geometry before initializing the temperature."
    @assert !isnan(thick3) "Layer 3 thickness (most likely thermal mantle lithosphere) is NaN. Make sure you initialized thermal mantle lithosphere thickness before initializing the temperature."

    term1 = 0.5 * hp1 / k1 * thick1 * thick1
    term2 = hp1 / k1 * thick1 * thick2 * k1 / k2 + 0.5*hp2 / k2 *thick2 * thick2
    term3 = (
        (hp1 / k1 * thick1 * thick3 * k1 / k2 + hp2 / k2 * thick2 * thick3) * k2 / k3 
        + 0.5 * hp3 / k3 * thick3 * thick3
        )
    integral1 = (term1 + term2 + term3) / thick_total
    k_harm    = 1.0/((thick1 / k1 + thick2 / k2 + thick3 / k3) / thick_total)

    local term_a, term_c

    if z_depth <= thick1
        term_a = 0.5 * hp1 / k1 * z_depth * z_depth
        term_c = z_depth/k1
    elseif thick1 + thick2 >= z_depth > thick1
        depth_term1 = z_depth - thick1
        term_a = term1 + hp1 / k1 * thick1 *depth_term1 *k1 / k2 + 0.5*hp2/k2*(depth_term1^2)
        term_c = thick1 / k1 + depth_term1 / k2
    elseif z_depth > thick1 + thick2
        depth_term2 = z_depth - thick1 - thick2
        term_a = (
            term1 + term2
            + (hp1 / k1 * thick1 * depth_term2 * k1 / k2 + hp2 /k2 * thick2 * depth_term2) * k2/ k3
            + 0.5 * hp3 / k3 * (depth_term2^2)
        )
        term_c = thick1 / k1 + thick2/k2 + (z_depth - thick1 - thick2) / k3
    end
    avg_grad = (temperature_bottom - temperature_top)/thick_total
    temperature = temperature_top - term_a + k_harm*(avg_grad + integral1)*term_c
    surface_heat_flow = k_harm*(avg_grad + integral1)/k1*k1*1000.0
    return temperature, surface_heat_flow
end

end # module 