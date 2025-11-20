module UpdateRockProps

using EarthBox.ModelDataContainer: ModelData
using EarthBox.GridFuncs: GridFuncs
using EarthBox.RockProperties.ThermalConductivityModel.UpdateManager.SekiguchiWaples: 
    thermal_conductivity_sekiguchi_waples

""" Calculate rock properties for all markers.

# Arguments
- `model::ModelData`: ModelData object
- `mantle_ids::Vector{Int16}`: Array of mantle material IDs
- `thermal_conductivity_ref_serpentinite::Float64`: Reference thermal conductivity for 
    serpentinite (W/m/K) used in Sekiguchi and Waples (2000) model
- `maximum_density_reduction::Float64`: Maximum density reduction (kg/m3) for 100% 
    serpentinization
- `characteristic_mantle_density::Float64`: Characteristic un-serpentinized mantle 
    density (kg/m3)

# Updated Arrays
- `model.markers.arrays.material.marker_rho`: Marker density in kg/m^3
- `model.markers.arrays.thermal.marker_rhocp`: Marker density multiplied by heat 
    capacity in J/m^3/K
- `model.markers.arrays.thermal.marker_kt`: Marker thermal conductivity in W/m/K
"""
function update_rock_props!(
    model::ModelData,
    mantle_ids::Vector{Int16},
    inside_flags::Vector{Int8};
    thermal_conductivity_ref_serpentinite::Float64 = 2.6,
    maximum_density_reduction::Float64 = 500.0,
    characteristic_mantle_density::Float64 = 3300.0
)::Nothing
    marker_TK = model.markers.arrays.thermal.marker_TK.array
    marker_kt = model.markers.arrays.thermal.marker_kt.array
    marker_rho = model.markers.arrays.material.marker_rho.array
    marker_rhocp = model.markers.arrays.thermal.marker_rhocp.array
    marker_matid = model.markers.arrays.material.marker_matid.array
    marker_serpentinization = model.markers.arrays.material.marker_serpentinization.array

    max_density_reduction_factor = calculate_maximum_serpentinization_density_reduction_factor(
        maximum_density_reduction, characteristic_mantle_density
    )

    marknum = model.markers.parameters.distribution.marknum.value
    for imarker in 1:marknum
        @inbounds matid = marker_matid[imarker]
        serpentinization_ratio = marker_serpentinization[imarker]
        if matid in mantle_ids && serpentinization_ratio > 0.0
            if inside_flags[imarker] == 1
                @inbounds begin
                    temperature_kelvins = marker_TK[imarker]
                    conductivity = marker_kt[imarker]
                    density = marker_rho[imarker]
                    rhocp = marker_rhocp[imarker]
                end

                serpentinite_conductivity = thermal_conductivity_sekiguchi_waples(
                    thermal_conductivity_ref_serpentinite,
                    temperature_kelvins
                )

                @inbounds serpentinization_ratio = marker_serpentinization[imarker]

                conductivity = (1.0 - serpentinization_ratio) * conductivity +
                             serpentinization_ratio * serpentinite_conductivity

                reduction_factor = calculate_density_reduction_factor(
                    serpentinization_ratio, max_density_reduction_factor)
                density = density * reduction_factor
                rhocp = rhocp * reduction_factor

                @inbounds begin
                    marker_kt[imarker] = conductivity
                    marker_rho[imarker] = density
                    marker_rhocp[imarker] = rhocp
                end
            end
        end
    end
    return nothing
end

""" Calculate maximum serpentinization density reduction factor.

# Arguments
- `maximum_density_reduction::Float64`: Maximum density reduction (kg/m3) for 100% 
    serpentinization
- `characteristic_mantle_density::Float64`: Characteristic un-serpentinized mantle 
    density (kg/m3)
"""
@inline function calculate_maximum_serpentinization_density_reduction_factor(
    maximum_density_reduction::Float64 = 500.0,
    characteristic_mantle_density::Float64 = 3300.0
)::Float64
    maximum_density_reduction_factor = (
        (characteristic_mantle_density - maximum_density_reduction) /
        characteristic_mantle_density
    )
    return maximum_density_reduction_factor
end

@inline function calculate_density_reduction_factor(
    serpentinization_ratio::Float64,
    maximum_density_reduction_factor::Float64
)::Float64
    reduction_factor = 1.0 + (maximum_density_reduction_factor - 1.0) * 
                      serpentinization_ratio
    return reduction_factor
end

end # module 