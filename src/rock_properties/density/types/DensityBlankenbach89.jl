module DensityBlankenbach89

using EarthBox.ModelDataContainer: ModelData

"""
    update_density!(model::ModelData)

Update density for all markers using the Blankenbach89 formulation.

# Arguments
- `model::ModelData`: The model data container containing marker information

# Updates
- `model.markers.arrays.material.marker_rho.array`: Array of marker densities
"""
function update_density!(model::ModelData, inside_flags::Vector{Int8})
    marker_rho = model.markers.arrays.material.marker_rho.array
    marker_matid = model.markers.arrays.material.marker_matid.array
    marknum = model.markers.parameters.distribution.marknum.value
    geometry = model.grids.parameters.geometry

    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                material_id = marker_matid[imarker]
                temperature_kelvins = model.markers.arrays.thermal.marker_TK.array[imarker]
                density_standard = model.materials.arrays.mat_rho.array[material_id, 1]
                expansivity = model.materials.arrays.mat_rho.array[material_id, 2]
                temperature_top_kelvins = model.bcs.parameters.temperature.temperature_top.value
            end
            density = marker_density_model(
                density_standard,
                expansivity,
                temperature_kelvins,
                temperature_top_kelvins
            )
            @inbounds marker_rho[imarker] = density
        end
    end
end

"""
    marker_density_model(
        density_standard::Float64,
        expansivity::Float64,
        temperature_kelvins::Float64,
        temperature_top_kelvins::Float64
    )::Float64

Calculate marker density using the Blankenbach89 formulation.

# Arguments
- `density_standard::Float64`: Standard density at reference conditions
- `expansivity::Float64`: Thermal expansivity coefficient
- `temperature_kelvins::Float64`: Current temperature in Kelvin
- `temperature_top_kelvins::Float64`: Reference temperature in Kelvin

# Returns
- `Float64`: Calculated density
"""
@inline function marker_density_model(
    density_standard::Float64,
    expansivity::Float64,
    temperature_kelvins::Float64,
    temperature_top_kelvins::Float64
)::Float64
    return density_blankenbach89(
        density_standard,
        expansivity,
        temperature_kelvins,
        temperature_top_kelvins
    )
end

"""
    density_blankenbach89(
        standard_density::Float64,
        thermal_expansion::Float64,
        temperature::Float64,
        temperature_ref::Float64
    )::Float64

Calculate density using formulation from Liao and Gerya (2014).

# Arguments
- `standard_density::Float64`: Density (kg/m^3) at standard reference temperature and pressure
- `thermal_expansion::Float64`: Thermal expansivity (1/K)
- `temperature::Float64`: Temperature (K)
- `temperature_ref::Float64`: Reference temperature (K)

# Returns
- `Float64`: Density (kg/m^3) at input pressure and temperature
"""
@inline function density_blankenbach89(
    standard_density::Float64,
    thermal_expansion::Float64,
    temperature::Float64,
    temperature_ref::Float64
)::Float64
    return standard_density * (1.0 - thermal_expansion * (temperature - temperature_ref))
end

end # module 