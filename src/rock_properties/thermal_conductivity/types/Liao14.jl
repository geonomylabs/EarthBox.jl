module Liao14

using EarthBox.ModelDataContainer: ModelData

"""
    update_conductivity!(model::ModelData)

Update thermal conductivity for all markers using the Liao14 formulation.

# Arguments
- `model::ModelData`: The model data container containing marker information

# Updates
- `model.markers.arrays.thermal.marker_kt.array`: Array of marker thermal conductivities
"""
function update_conductivity!(model::ModelData, inside_flags::Vector{Int8})
    marker_kt = model.markers.arrays.thermal.marker_kt.array
    marker_matid = model.markers.arrays.material.marker_matid.array
    marknum = model.markers.parameters.distribution.marknum.value

    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                material_id = marker_matid[imarker]
                pressure_pascals = model.markers.arrays.pressure.marker_pr.array[imarker]
                temperature_kelvins = model.markers.arrays.thermal.marker_TK.array[imarker]
                conductivity_reference = model.materials.arrays.mat_kt.array[material_id, 1]
                conductivity_factor_a = model.materials.arrays.mat_kt.array[material_id, 2]
            end
            conductivity = thermal_conductivity_liao14(
                conductivity_reference,
                conductivity_factor_a,
                pressure_pascals,
                temperature_kelvins
            )
            @inbounds marker_kt[imarker] = conductivity
        end
    end
end

""" Calculate thermal conductivity using the Liao14 formulation.

Formulation from Liao and Gerya (2014):

    k = [k0 + a/(T + 77)]*exp(4e-5*P)

where k0 is the reference thermal conductivity (W/m/K), a is a term in the
thermal conductivity equation, T is the temperature (K), and P is the
pressure (MPa). This formulation is based on the work of Clauser and
Huenges (1995) and Zoth and Hanel (1988). Note that the pressure term was
added by Liao and Gerya (2014) based on data from Clauser and Huenges
(1995). Note that this equation is modified for temperature in Kelvin
whereas the original formulation uses Celsius.

# Arguments
- `thermal_conductivity_k0::Float64`: Reference thermal conductivity (W/m/K)
- `thermal_conductivity_a::Float64`: Term a in thermal conductivity equation
- `pressure::Float64`: Pressure (Pa)
- `temperature::Float64`: Temperature (K)

# Returns
- `Float64`: Thermal conductivity (W/m/K) at input pressure and temperature
"""
@inline function thermal_conductivity_liao14(
    thermal_conductivity_k0::Float64,
    thermal_conductivity_a::Float64,
    pressure::Float64,
    temperature::Float64
)::Float64
    return (thermal_conductivity_k0 + thermal_conductivity_a / (temperature + 77.0)) *
           exp(0.00004 * pressure / 1e6)
end

end # module 