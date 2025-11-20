module Uniform

using EarthBox.ModelDataContainer: ModelData

"""
    update_conductivity!(model::ModelData)

Update thermal conductivity for all markers using a uniform conductivity model.

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
                conductivity_reference = model.materials.arrays.mat_kt.array[material_id, 1]
                marker_kt[imarker] = conductivity_reference
            end
        end
    end
end

end # module 