module TemperatureDependent

using EarthBox.ModelDataContainer: ModelData
using EarthBox.GridFuncs: check_in_domain

"""
    update_conductivity!(model::ModelData)

Update thermal conductivity for all markers using a temperature-dependent conductivity model.

# Arguments
- `model::ModelData`: The model data container containing marker information

# Updates
- `model.markers.arrays.thermal.marker_kt.array`: Array of marker thermal conductivities
"""
function update_conductivity!(model::ModelData, inside_flags::Vector{Int8})
    marker_kt = model.markers.arrays.thermal.marker_kt.array
    marker_matid = model.markers.arrays.material.marker_matid.array
    marknum = model.markers.parameters.distribution.marknum.value
    temperature_top_kelvins = model.bcs.parameters.temperature.temperature_top.value

    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                material_id = marker_matid[imarker]
                temperature_kelvins = model.markers.arrays.thermal.marker_TK.array[imarker]
                conductivity_reference = model.materials.arrays.mat_kt.array[material_id, 1]
                conductivity_factor_a = model.materials.arrays.mat_kt.array[material_id, 2]
            end
            conductivity = marker_conductivity_model(
                conductivity_reference,
                conductivity_factor_a,
                temperature_kelvins,
                temperature_top_kelvins
            )
            @inbounds marker_kt[imarker] = conductivity
        end
    end
end

"""
    marker_conductivity_model(
        conductivity_reference::Float64,
        conductivity_factor_a::Float64,
        temperature_kelvins::Float64,
        temperature_top_kelvins::Float64
    )::Float64

Calculate marker thermal conductivity using the temperature-dependent model.

# Arguments
- `conductivity_reference::Float64`: Reference thermal conductivity
- `conductivity_factor_a::Float64`: Term a in thermal conductivity equation
- `temperature_kelvins::Float64`: Temperature (K)
- `temperature_top_kelvins::Float64`: Temperature at the top of model (K)

# Returns
- `Float64`: Calculated thermal conductivity
"""
@inline function marker_conductivity_model(
    conductivity_reference::Float64,
    conductivity_factor_a::Float64,
    temperature_kelvins::Float64,
    temperature_top_kelvins::Float64
)::Float64
    return thermal_conductivity_channel_flow_benchmark(
        conductivity_reference,
        conductivity_factor_a,
        temperature_kelvins,
        temperature_top_kelvins
    )
end

"""
    thermal_conductivity_channel_flow_benchmark(
        thermal_conductivity_k0::Float64,
        thermal_conductivity_a::Float64,
        temperature::Float64,
        temperature_top::Float64
    )::Float64

Calculate thermal conductivity for channel flow benchmark.

# Arguments
- `thermal_conductivity_k0::Float64`: Reference thermal conductivity (W/m/K)
- `thermal_conductivity_a::Float64`: Term a in thermal conductivity equation
- `temperature::Float64`: Temperature (K)
- `temperature_top::Float64`: Temperature at the top of model (K)

# Returns
- `Float64`: Thermal conductivity (W/m/K)
"""
@inline function thermal_conductivity_channel_flow_benchmark(
    thermal_conductivity_k0::Float64,
    thermal_conductivity_a::Float64,
    temperature::Float64,
    temperature_top::Float64
)::Float64
    return thermal_conductivity_k0 / (
        1.0 + thermal_conductivity_a * (temperature - temperature_top) / temperature_top
    )
end

end # module 