module Blankenback89

import EarthBox.ModelDataContainer: ModelData
import EarthBox.GridFuncs: check_in_domain
import Base.Threads
import ..FlowUtils: apply_minimum_temperature_limit

function update_marker_flow_viscosity!(model::ModelData, inside_flags::Vector{Int8})::Nothing
    temperature_top = model.bcs.parameters.temperature.temperature_top.value
    temperature_bottom = model.bcs.parameters.temperature.temperature_bottom.value
    ysize = model.grids.parameters.geometry.ysize.value
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marker_TK = model.markers.arrays.thermal.marker_TK.array
    marker_eta_flow = model.markers.arrays.rheology.marker_eta_flow.array
    marker_preexp = model.markers.arrays.rheology.marker_preexp.array
    mat_flow = model.materials.arrays.mat_flow.array
    marker_matid = model.markers.arrays.material.marker_matid.array
    geometry = model.grids.parameters.geometry
    marknum = model.markers.parameters.distribution.marknum.value

    Threads.@threads for imarker in 1:marknum
        x_marker = marker_x[imarker]
        y_marker = marker_y[imarker]
        if check_in_domain(geometry, x_marker, y_marker)
            mat_id = marker_matid[imarker]
            temperature = apply_minimum_temperature_limit(
                marker_TK[imarker], temperature_top)
            viscosity_ref = mat_flow[mat_id,15]
            b_term = mat_flow[mat_id,16]
            c_term = mat_flow[mat_id,17]
            viscosity = calculate_viscosity_blankenbach89(
                y_marker,
                ysize,
                temperature_top,
                temperature_bottom,
                viscosity_ref,
                b_term,
                c_term,
                temperature
            )
            marker_eta_flow[imarker] = viscosity
            marker_preexp[imarker] = viscosity_ref
        end
    end
    return nothing
end

function calculate_viscosity_blankenbach89(
    y_location::Float64,
    ysize::Float64,
    temperature_top::Float64,
    temperature_bottom::Float64,
    viscosity_ref::Float64,
    b_term::Float64,
    c_term::Float64,
    temperature::Float64
)::Float64
    viscosity = viscosity_ref * exp(
        -b_term * (temperature - temperature_top) /
        (temperature_bottom - temperature_top)
        + c_term * y_location / ysize
    )
    return viscosity
end

end # module
