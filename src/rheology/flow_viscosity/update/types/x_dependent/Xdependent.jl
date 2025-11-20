module Xdependent

import EarthBox.ModelDataContainer: ModelData
import EarthBox.GridFuncs: check_in_domain
import Base.Threads
import ..FlowUtils: apply_minimum_temperature_limit

"""
    update_marker_flow_viscosity!(model::ModelData)::Nothing

Update marker viscosity and pre-exponential arrays for x-dependent flow law.
"""
function update_marker_flow_viscosity!(model::ModelData, inside_flags::Vector{Int8})::Nothing
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marker_matid = model.markers.arrays.material.marker_matid.array
    mat_flow = model.materials.arrays.mat_flow.array
    marker_eta_flow = model.markers.arrays.rheology.marker_eta_flow.array
    marker_preexp = model.markers.arrays.rheology.marker_preexp.array
    viscosity_maximum = model.materials.parameters.viscosity_limits.viscosity_max.value
    geometry = model.grids.parameters.geometry
    marknum = model.markers.parameters.distribution.marknum.value

    Threads.@threads for imarker in 1:marknum
        x_marker = marker_x[imarker]
        y_marker = marker_y[imarker]
        if check_in_domain(geometry, x_marker, y_marker)
            mat_id = marker_matid[imarker]
            viscosity = FlowLaws.calculate_viscosity_x_dependent(x_marker, viscosity_maximum)
            marker_eta_flow[imarker] = viscosity
            marker_preexp[imarker] = mat_flow[mat_id][1]
        end
    end
    return nothing
end

"""
    calculate_viscosity_x_dependent(x_location::Float64, viscosity_maximum::Float64)::Float64

Simple viscosity model with lateral variation.
"""
function calculate_viscosity_x_dependent(
    x_location::Float64,
    viscosity_maximum::Float64
)::Float64
    pressure_top = 1e9
    xdim = 10000.0
    ydim = 9500.0

    n_value = 3.0
    c_value = 1e-37

    pressure_gradient = -pressure_top / ydim

    if abs(x_location - xdim/2.0) > 0.0
        viscosity = (
            1.0/c_value
            * (-pressure_gradient)^(1-n_value)
            * (x_location - xdim/2.0)^(1 - n_value)
        )
    else
        viscosity = viscosity_maximum
    end
    return viscosity
end

end # module

