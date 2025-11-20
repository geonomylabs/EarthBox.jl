module TemperatureDependentWaples

using EarthBox.ModelDataContainer: ModelData
using EarthBox.ConversionFuncs: kelvin_to_celsius

"""
    update_rhocp!(model::ModelData)

Calculate rhocp for all markers using variable heat capacity based on temperature.

# Arguments
- `model::ModelData`: The model data container containing marker information

# Updates
- `model.markers.arrays.thermal.marker_rhocp.array`: Array of marker rhocp values
"""
function update_rhocp!(model::ModelData, inside_flags::Vector{Int8})
    marker_rho = model.markers.arrays.material.marker_rho.array
    marker_rhocp = model.markers.arrays.thermal.marker_rhocp.array
    marker_matid = model.markers.arrays.material.marker_matid.array
    mat_cp = model.materials.arrays.mat_cp.array
    marknum = model.markers.parameters.distribution.marknum.value
    marker_TK = model.markers.arrays.thermal.marker_TK.array
    maximum_heat_capacity = model.heat_equation.parameters.rhocp.maximum_heat_capacity.value

    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                material_id = marker_matid[imarker]
                temperature_celsius = kelvin_to_celsius(marker_TK[imarker])
                heat_capacity_at_20celsius = mat_cp[material_id]
                heat_capacity = calculate_heat_capacity_waples(
                    temperature_celsius,
                    heat_capacity_at_20celsius
                )
                heat_capacity = min(heat_capacity, maximum_heat_capacity)
                density = marker_rho[imarker]
            end
            rhocp = density * heat_capacity
            @inbounds marker_rhocp[imarker] = rhocp
        end
    end
end

"""
    calculate_heat_capacity_waples(
        temperature_celsius::Float64,
        heat_capacity_at_20celsius::Float64
    )::Float64

Calculate heat capacity using Waples (1985) model.

See Hantschel and Kauerauf (2010) for details.

A cut off of 1200 J/kg/K is applied to avoid heat capacity that avoid
problems associated with this empirical law at high temperatures.

This issue was encountered when modeling melting at mid ocean ridges. If the
cut off was not included, the model would produce heat capacity ranging from
1250 to 1400 J/kg/K in the upper mantle yielding higher temperature and too
much melt (and too much crustal thickness) at normal spreading centers with 
bottom temperatures of 1330C.

# Arguments
- `temperature_celsius::Float64`: Temperature in Celsius
- `heat_capacity_at_20celsius::Float64`: Heat capacity at 20 Celsius (J/kg/K)

# Returns
- `Float64`: Heat capacity in J/kg/K
"""
@inline function calculate_heat_capacity_waples(
    temperature_celsius::Float64,
    heat_capacity_at_20celsius::Float64
)::Float64
    return heat_capacity_at_20celsius * (
        0.953 +
        2.29e-3 * temperature_celsius -
        2.835e-6 * temperature_celsius^2 +
        1.191e-9 * temperature_celsius^3
    )
end

end # module 