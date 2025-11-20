module ConstantRhoCp

using EarthBox.ModelDataContainer: ModelData

"""
    update_rhocp!(model::ModelData)

Calculate rhocp for all markers using constant heat capacity.

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

    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                material_id = marker_matid[imarker]
                heat_capacity = mat_cp[material_id]
                density = marker_rho[imarker]
                rhocp = density * heat_capacity
            end
            @inbounds marker_rhocp[imarker] = rhocp
        end
    end
end

end # module 