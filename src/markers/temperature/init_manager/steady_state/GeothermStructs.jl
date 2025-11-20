module GeothermStructs

"""
Struct for the steady-state analytical 3 layer geotherm.
"""
struct SteadyState3LayerData
    # Temperature at the sticky-rock interface of the model in Kelvin
    temperature_top::Float64
    # Temperature at the bottom of the model in Kelvin
    temperature_bottom::Float64
    # Thickness of the first layer in meters (e.g. upper crust)
    layer1_thickness::Float64
    # Thickness of the second layer in meters (e.g. lower crust)
    layer2_thickness::Float64
    # Thickness of the third layer in meters (e.g. thermal mantle lithosphere)
    layer3_thickness::Float64
    # Total thickness of the thermal mantle lithosphere in meters
    thickness_total::Float64
    # Heat production of the first layer in W/m^3
    layer1_heat_production::Float64
    # Heat production of the second layer in W/m^3
    layer2_heat_production::Float64
    # Heat production of the third layer in W/m^3
    layer3_heat_production::Float64
    # Thermal conductivity of the first layer in W/mK
    layer1_conductivity::Float64
    # Thermal conductivity of the second layer in W/mK
    layer2_conductivity::Float64
    # Thermal conductivity of the third layer in W/m/K
    layer3_conductivity::Float64
end

end # module 