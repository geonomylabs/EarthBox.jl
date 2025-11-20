module DensityLiao14

using EarthBox.ModelDataContainer: ModelData
using EarthBox.GridFuncs: check_in_domain

""" Update density for all markers using the Liao14 model.

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
                pressure_pascals = model.markers.arrays.pressure.marker_pr.array[imarker]
                temperature_kelvins = model.markers.arrays.thermal.marker_TK.array[imarker]
                density_standard = model.materials.arrays.mat_rho.array[material_id, 1]
                expansivity = model.materials.arrays.mat_rho.array[material_id, 2]
                compressibility = model.materials.arrays.mat_rho.array[material_id, 3]
            end
            
            density = marker_density_model(
                density_standard,
                expansivity,
                compressibility,
                pressure_pascals,
                temperature_kelvins
            )
            @inbounds marker_rho[imarker] = density
        end
    end
end

""" Calculate marker density using reference conditions.

# Arguments
- `density_standard`: Standard density (kg/m³)
- `expansivity`: Thermal expansivity (1/K)
- `compressibility`: Compressibility (1/Pa)
- `pressure_pascals`: Pressure (Pa)
- `temperature_kelvins`: Temperature (K)

# Returns
- `Float64`: Calculated density (kg/m³)
"""
@inline function marker_density_model(
    density_standard::Float64,
    expansivity::Float64,
    compressibility::Float64,
    pressure_pascals::Float64,
    temperature_kelvins::Float64
)::Float64
    temperature_ref = 298.15
    pressure_ref = 1e5
    density = density_liao14(
        density_standard,
        expansivity,
        compressibility,
        pressure_pascals,
        temperature_kelvins,
        temperature_ref,
        pressure_ref
    )
    return density
end

""" Calculate density using formulation from Liao and Gerya (2014).

# Arguments
- `standard_density`: Density (kg/m³) at standard reference temperature and pressure
- `thermal_expansion`: Thermal expansivity (1/K)
- `compressibility`: Compressibility (1/Pa)
- `pressure`: Pressure (Pa)
- `temperature`: Temperature (K)
- `temperature_ref`: Reference temperature (K)
- `pressure_ref`: Reference pressure (Pa)

# Returns
- `Float64`: Density (kg/m³) at input pressure and temperature
"""
@inline function density_liao14(
    standard_density::Float64,
    thermal_expansion::Float64,
    compressibility::Float64,
    pressure::Float64,
    temperature::Float64,
    temperature_ref::Float64,
    pressure_ref::Float64
)::Float64
    density = standard_density *
        (1.0 - thermal_expansion * (temperature - temperature_ref)) *
        (1.0 + compressibility * (pressure - pressure_ref))
    return density
end

end # module DensityLiao14 