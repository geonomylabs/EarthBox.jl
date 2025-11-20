module TransientBottomCalculator

""" Calculate thermal boundary conditions for transient cooling.

Transient cooling involves an initial temperature perturbation to the
base of the lithosphere and bottom of model at time t=0 (warmer state if 
delta is +). The perturbation at the base of the model is turned off during 
a user specified time interval.

This is used to simulate an initially warm asthenosphere associated with
a plume that spreads out at the base of the lithosphere and cools over time.

# Arguments
- `temperature_base_lith::Float64`: Cooler transient temperature at the base of the 
    lithosphere in kelvin
- `delta_temperature_transient::Float64`: Initial temperature perturbation in delta 
    kelvin or delta celsius
- `adiabatic_gradient::Float64`: Adiabatic temperature gradient in celsius/km or kelvin/km
- `ysize::Float64`: Size of model in y-direction in meters
- `thick_air::Float64`: Thickness of sticky air layer in meters
- `thick_lith::Float64`: Thickness of lithosphere in meters

# Returns
- `temperature_base_lith_warmer::Float64`: Initial elevated temperature 
  at the base of the lithosphere in Kelvin.
- `temperature_bottom_warmer::Float64`: Initial elevated temperature 
  at the bottom of the model in Kelvin
- `temperature_bottom_cooler::Float64`: Cooler transient temperature 
  at the bottom of the model in Kelvin

# Notes
- 'temperature_base_lith' is should be compared to the AnalyticalThreeLayer model
for consistency. Currently this parameter is being used as a reference point for
defining the bottom temperature. 
"""
function calculate_thermal_bcs_for_transient_cooling(
    temperature_base_lith::Float64,
    delta_temperature_transient::Float64,
    adiabatic_gradient::Float64,
    ysize::Float64,
    thick_air::Float64,
    thick_lith::Float64
)::Tuple{Float64, Float64, Float64}
    distance_from_base_lith_to_bottom = 
        calculate_distance_from_base_lithosphere_to_bottom(ysize, thick_air, thick_lith)
    temperature_base_lith_warmer = temperature_base_lith + delta_temperature_transient
    temperature_bottom_warmer = calculate_bottom_temperature(
        temperature_base_lith_warmer,
        adiabatic_gradient,
        distance_from_base_lith_to_bottom
		)
    temperature_bottom_cooler = calculate_bottom_temperature(
        temperature_base_lith,
        adiabatic_gradient,
        distance_from_base_lith_to_bottom
		)
    return (
		temperature_base_lith_warmer,
        temperature_bottom_warmer,
        temperature_bottom_cooler
		)
end

""" Calculate distance from base of lithosphere to bottom of model.

Assumes length unit are consistent for all input variables.

# Arguments
- `ysize::Float64`: Size of model in y-direction
- `thick_air::Float64`: Thickness of sticky air layer
- `thick_lithosphere::Float64`: Thickness of lithosphere

# Returns
- `distance_from_base_lithosphere_to_bottom::Float64`: Distance from base of 
  lithosphere to bottom of model in input variable units
"""
function calculate_distance_from_base_lithosphere_to_bottom(
    ysize::Float64,
    thick_air::Float64,
    thick_lithosphere::Float64
)::Float64
    return ysize - thick_air - thick_lithosphere
end

""" Calculate base lithosphere temperature taking into account transient 
temperature perturbation.
"""
function calculate_temperature_base_lithosphere(
        temperature_base_lith::Float64,
        delta_temperature_transient::Float64=0.0
)::Float64
    temperature_base_lith += delta_temperature_transient
    return temperature_base_lith
end

function calculate_bottom_temperature(
        temperature_base_lith::Float64,
        adiabatic_gradient::Float64,
        distance_from_base_lith_to_bottom::Float64,
        delta_temperature_transient::Float64=0.0
)::Float64
    temperature_bottom = (
        temperature_base_lith
        + delta_temperature_transient
        + distance_from_base_lith_to_bottom
        * adiabatic_gradient / 1000.0)
    return temperature_bottom
end
end