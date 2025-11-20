module TempIniStructs

"""
Immutable data structures for initial temperature condition calculations.
"""

"""
Struct for the thermal conductivity of the three layers in the analytical
three layer model.
"""
struct Conductivity
    conductivity_upper_crust::Float64  # W/mK
    conductivity_lower_crust::Float64  # W/mK
    conductivity_mantle::Float64       # W/mK
end

"""
Struct for heat production used in the steady-state analytical 3 layer
geotherm model.
"""
struct HeatProduction
    heat_production_upper_crust::Float64  # W/m^3
    heat_production_lower_crust::Float64  # W/m^3
    heat_production_mantle::Float64       # W/m^3
end

"""
Struct for the perturbation used in the analytical three layer model
and the four linear segments model.
"""
struct Perturbation
    amplitude_perturbation::Float64  # meters
    width_perturbation::Float64      # meters
end

"""
Temperature boundary condition data structure used in both the analytical
three layer model and the four linear segments model.
"""
struct TemperatureBCs
    temperature_top::Float64        # Temperature at the top of the model in Kelvin
    temperature_bottom::Float64     # Temperature at the base of the model in Kelvin
    temperature_base_lith::Float64  # Temperature at the base of the lithosphere in Kelvin
end

"""
Struct for the four linear segments model.
"""
struct InternalTemperature
    temperature_surface::Float64  # Temperature at the sticky-rock interface in Kelvin
    temperature_moho::Float64     # Temperature at the Moho in Kelvin
end

"""
Struct for layer thickness information.
"""
struct LayerThickness
    thick_air::Float64              # meters
    thick_upper_crust::Float64      # meters
    thick_lower_crust::Float64      # meters
    thick_upper_lith::Float64       # meters
    thick_middle_lith::Float64      # meters
    thick_lower_lith::Float64       # meters
    thick_thermal_lithosphere::Float64  # meters
end

export Conductivity, HeatProduction, Perturbation, TemperatureBCs, InternalTemperature, LayerThickness

end # module