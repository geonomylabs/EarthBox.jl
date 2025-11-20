module DensityNaliboff17

using EarthBox.ModelDataContainer: ModelData
using EarthBox.GridFuncs: check_in_domain

"""
    update_density!(model::ModelData)

Update density for all markers using the Naliboff17 model.

# Arguments
- `model::ModelData`: The model data container containing marker information

# Updates
- `model.markers.arrays.material.marker_rho.array`: Array of marker densities
"""
function update_density!(model::ModelData, inside_flags::Vector{Int8})
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
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
            end
            
            density = density_naliboff17(
                density_standard,
                expansivity,
                temperature_kelvins
            )
            @inbounds marker_rho[imarker] = density
        end
    end
end

"""
    density_naliboff17(
        standard_density::Float64,
        thermal_expansivity::Float64,
        temperature::Float64,
        temperature_ref::Float64 = 830.0
    )::Float64

Calculate density using formulation from Naliboff et al. 2017.

# Arguments
- `standard_density`: Density (kg/m³) at standard reference temperature
- `thermal_expansivity`: Thermal expansivity (1/K)
- `temperature`: Temperature (K)
- `temperature_ref`: Reference temperature (K), defaults to 830.0

# Returns
- `Float64`: Density (kg/m³) at input temperature
"""
@inline function density_naliboff17(
    standard_density::Float64,
    thermal_expansivity::Float64,
    temperature::Float64,
    temperature_ref::Float64 = 830.0
)::Float64
    density = standard_density * (1.0 - thermal_expansivity * (temperature - temperature_ref))
    return density
end

end # module DensityNaliboff17 