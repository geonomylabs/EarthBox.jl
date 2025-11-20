module CouetteFlowBenchmark

import EarthBox.ModelDataContainer: ModelData
import EarthBox.GridFuncs: check_in_domain
import Base.Threads
import ..FlowUtils: apply_minimum_temperature_limit

function update_marker_flow_viscosity!(model::ModelData, inside_flags::Vector{Int8})::Nothing
    temperature_left = model.bcs.parameters.temperature.temperature_left.value
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marker_TK = model.markers.arrays.thermal.marker_TK.array
    marker_eta_flow = model.markers.arrays.rheology.marker_eta_flow.array
    marker_preexp = model.markers.arrays.rheology.marker_preexp.array
    mat_flow = model.materials.arrays.mat_flow.array
    marker_matid = model.markers.arrays.material.marker_matid.array
    gas_constant = model.RGAS.value
    geometry = model.grids.parameters.geometry
    marknum = model.markers.parameters.distribution.marknum.value

    Threads.@threads for imarker in 1:marknum
        x_marker = marker_x[imarker]
        y_marker = marker_y[imarker]
        if check_in_domain(geometry, x_marker, y_marker)
            mat_id = marker_matid[imarker]
            temperature = apply_minimum_temperature_limit(marker_TK[imarker], temperature_left)
            pre_exponential = mat_flow[mat_id, 13]
            activation_energy = mat_flow[mat_id, 14]
            temperature_left_wall = temperature_left
            viscosity = calculate_viscosity_couette_flow_benchmark(
                temperature_left_wall,
                pre_exponential,
                activation_energy,
                gas_constant,
                temperature
            )
            marker_eta_flow[imarker] = viscosity
            marker_preexp[imarker] = pre_exponential
        end
    end
    return nothing
end

function calculate_viscosity_couette_flow_benchmark(
    temperature_left_wall::Float64,
    pre_exponential::Float64,
    activation_energy::Float64,
    gas_constant::Float64,
    temperature::Float64
)::Float64
    viscosity = pre_exponential * exp(
        activation_energy * 1000.0 / gas_constant / temperature_left_wall
        * (1.0 - (temperature - temperature_left_wall) / temperature_left_wall)
    )
    return viscosity
end

end # module

