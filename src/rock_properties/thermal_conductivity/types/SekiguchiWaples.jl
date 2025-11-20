module SekiguchiWaples

using EarthBox.ModelDataContainer: ModelData

"""
    update_conductivity!(model::ModelData)

Update thermal conductivity for all markers using the Sekiguchi-Waples formulation.

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
                temperature_kelvins = model.markers.arrays.thermal.marker_TK.array[imarker]
                conductivity_reference = model.materials.arrays.mat_kt.array[material_id, 1]
            end
            
            # Only use the Sekiguchi-Waples model if the reference
            # conductivity is less than or equal to 20.0 W/m/K. Otherwise,
            # use the reference conductivity since these cases involve
            # geologically unrealistic thermal conductivities used to mimic
            # a sediment-water interface.
            if conductivity_reference <= 20.0
                conductivity = thermal_conductivity_sekiguchi_waples(
                    conductivity_reference,
                    temperature_kelvins
                )
            else
                conductivity = conductivity_reference
            end
            @inbounds marker_kt[imarker] = conductivity
        end
    end
end

"""
    thermal_conductivity_sekiguchi_waples(
        thermal_conductivity_k0::Float64,
        temperature::Float64
    )::Float64

Calculate thermal conductivity using the Sekiguchi-Waples model.

Matrix thermal conductivity model adapted from Sekiguchi (1984) and
modified by Doug Waples as summarized by Hantschel and Kauerauf (2009).

# Arguments
- `thermal_conductivity_k0::Float64`: Reference thermal conductivity at 20 degrees Celsius (W/m/K)
- `temperature::Float64`: Temperature (K)

# Returns
- `Float64`: Matrix thermal conductivity (W/m/K) at input temperature
"""
@inline function thermal_conductivity_sekiguchi_waples(
    thermal_conductivity_k0::Float64,
    temperature::Float64
)::Float64
    return 358.0 * (1.0227 * thermal_conductivity_k0 - 1.882) * 
           (1.0 / temperature - 0.00068) + 1.84
end

end # module 