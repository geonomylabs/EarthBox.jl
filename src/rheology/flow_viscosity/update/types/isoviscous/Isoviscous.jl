module Isoviscous

import EarthBox.ModelDataContainer: ModelData
import EarthBox.GridFuncs: check_in_domain
import EarthBox.Arrays: ArrayUtils
import Base.Threads
import ..FlowUtils: apply_minimum_temperature_limit

function update_marker_flow_viscosity!(model::ModelData, inside_flags::Vector{Int8})::Nothing
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marker_matid = model.markers.arrays.material.marker_matid.array
    marker_eta_flow = model.markers.arrays.rheology.marker_eta_flow.array
    marker_preexp = model.markers.arrays.rheology.marker_preexp.array
    mat_flow = model.materials.arrays.mat_flow.array
    geometry = model.grids.parameters.geometry
    marknum = model.markers.parameters.distribution.marknum.value

    Threads.@threads for imarker in 1:marknum
        x_marker = marker_x[imarker]
        y_marker = marker_y[imarker]
        if check_in_domain(geometry, x_marker, y_marker)
            mat_id = marker_matid[imarker]
            viscosity = mat_flow[mat_id, 1]
            marker_eta_flow[imarker] = viscosity
            marker_preexp[imarker] = viscosity
        end
    end

    return nothing
end

end # module

