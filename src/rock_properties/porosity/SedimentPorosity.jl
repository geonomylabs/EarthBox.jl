module SedimentPorosity

import EarthBox.ModelDataContainer: ModelData
import EarthBox: SedimentWaterInterface

""" Calculate sediment properties taking compaction into account.

# Arguments
- `heat_capacity`: Heat capacity of sediment in J/K/kg
- `conductivity`: Thermal conductivity of sediment in W/m/K
- `density`: Density of sediment in kg/m^3
- `porosity_at_mudline`: Porosity at mudline in fraction
- `porosity_decay_depth`: Decay depth of porosity in meters
- `conductivity_water`: Thermal conductivity of water in W/m/K
- `density_water`: Density of water in kg/m^3
- `heat_capacity_water`: Heat capacity of water in J/K/kg
- `x_location`: X location of marker in meters
- `y_location`: Y location of marker in meters
- `gridt`: Topography grid array with Y-coordinates of topography stored in grid[:,2]
- `max_burial`: Maximum burial depth of marker in meters

# Returns
- `conductivity`: Bulk thermal conductivity of sediment in W/M/K
- `density`: Bulk density of sediment in kg/m^3
- `rhocp`: Bulk density and heat capacity of sediment
"""
@inline function calculate_rock_props(
    heat_capacity::Float64,
    conductivity::Float64,
    density::Float64,
    porosity_at_mudline::Float64,
    porosity_decay_depth::Float64,
    conductivity_water::Float64,
    density_water::Float64,
    heat_capacity_water::Float64,
    x_location::Float64,
    y_location::Float64,
    gridt::Array{Float64,2},
    max_burial::Float64
)::Tuple{Float64,Float64,Float64}
    porosity = calculate_porosity(
        porosity_at_mudline, porosity_decay_depth,
        x_location, y_location, gridt, max_burial
    )
    bulk_conductivity = calculate_bulk_property(
        conductivity, conductivity_water, porosity)
    bulk_density = calculate_bulk_property(
        density, density_water, porosity)
    bulk_heat_capacity = calculate_bulk_property(
        heat_capacity, heat_capacity_water, porosity)
    bulk_rhocp = bulk_density * bulk_heat_capacity
    return bulk_conductivity, bulk_density, bulk_rhocp
end

""" Calculate porosity of sediment.

# Arguments
- `porosity_at_mudline`: Porosity at mudline in fraction
- `porosity_decay_depth`: Decay depth of porosity in meters
- `x_location`: X-location of marker in meters
- `y_location`: Y-location of marker in meters
- `gridt`: Topography grid array
- `max_burial`: Maximum burial depth of marker in meters

# Returns
- `porosity`: Porosity of sediment in fraction
"""
@inline function calculate_porosity(
    porosity_at_mudline::Float64,
    porosity_decay_depth::Float64,
    x_location::Float64,
    y_location::Float64,
    gridt::Array{Float64,2},
    max_burial::Float64
)::Float64
    y_mudline = SedimentWaterInterface.get_depth(x_location, gridt)
    y_submud = y_location - y_mudline
    y_submud = max(y_submud, max_burial)
    depth_decay_term = 1.0 / porosity_decay_depth
    porosity = porosity_at_mudline * exp(-depth_decay_term * y_submud)
    return porosity
end

@inline function calculate_bulk_property(
    matrix_property::Float64,
    water_property::Float64,
    porosity::Float64
)::Float64
    return matrix_property * (1.0 - porosity) + water_property * porosity
end 

end