module Naliboff17

using EarthBox.ModelDataContainer: ModelData

"""
    update_conductivity!(model::ModelData)

Update thermal conductivity for all markers using the Naliboff17 formulation.

# Arguments
- `model::ModelData`: The model data container containing marker information

# Updated arrays from group `model.markers.arrays.thermal`
- `marker_kt.array::Vector{Float64}`: Array of marker thermal conductivities
"""
function update_conductivity!(model::ModelData, inside_flags::Vector{Int8})
    marker_kt = model.markers.arrays.thermal.marker_kt.array
    marker_matid = model.markers.arrays.material.marker_matid.array
    marknum = model.markers.parameters.distribution.marknum.value

    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                material_id = marker_matid[imarker]
                y_location = model.markers.arrays.location.marker_y.array[imarker]
                conductivity_reference = model.materials.arrays.mat_kt.array[material_id, 1]
                conductivity = marker_conductivity_model(conductivity_reference, y_location)
            end
            @inbounds marker_kt[imarker] = conductivity
        end
    end
end

"""
    marker_conductivity_model(
        conductivity_reference::Float64,
        y_location::Float64
    )::Float64

Calculate marker thermal conductivity using the Naliboff17 formulation.

# Arguments
- `conductivity_reference::Float64`: Reference thermal conductivity
- `y_location::Float64`: Depth of marker in meters

# Returns
- `Float64`: Calculated thermal conductivity
"""
@inline function marker_conductivity_model(
    conductivity_reference::Float64,
    y_location::Float64
)::Float64
    return thermal_conductivity_naliboff17(y_location, conductivity_reference)
end

"""
    thermal_conductivity_naliboff17(
        y_location::Float64,
        thermal_conductivity::Float64;
        lithosphere_asthenosphere_boundary_depth::Float64=110000.0
    )::Float64

Increase thermal conductivity above 1300C to approximate an adiabatic gradient.

# Arguments
- `y_location::Float64`: Depth of marker in meters
- `thermal_conductivity::Float64`: Thermal conductivity (W/m/K)
- `lithosphere_asthenosphere_boundary_depth::Float64`: Depth of lithosphere-asthenosphere boundary in meters

# Returns
- `Float64`: Thermal conductivity (W/m/K)
"""
@inline function thermal_conductivity_naliboff17(
    y_location::Float64,
    thermal_conductivity::Float64;
    lithosphere_asthenosphere_boundary_depth::Float64=110000.0
)::Float64
    if y_location > lithosphere_asthenosphere_boundary_depth
        thermal_conductivity = 39.25
    end
    return thermal_conductivity
end

end # module 