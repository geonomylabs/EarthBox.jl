module Channel

import EarthBox.ModelDataContainer: ModelData

function recycle_markers!(model::ModelData)::Nothing
    execute_recycle_steps(model)
    return nothing
end

"""
    execute_recycle_steps(model::ModelData)::Nothing

Recycle markers for channel models.

# Updated Arrays
- `model.markers.arrays.location.marker_y.array`: y-location of markers (m)
- `model.markers.arrays.thermal.marker_TK.array`: Temperature of markers (K)
"""
function execute_recycle_steps(model::ModelData)::Nothing
    t_top = model.bcs.parameters.temperature.temperature_top.value
    t_bottom = model.bcs.parameters.temperature.temperature_bottom.value
    ynum = model.grids.parameters.geometry.ynum.value

    marker_y = model.markers.arrays.location.marker_y.array
    marker_temperature = model.markers.arrays.thermal.marker_TK.array

    gridy_min = model.grids.arrays.basic.gridy_b.array[1]
    gridy_max = model.grids.arrays.basic.gridy_b.array[ynum]
    ysize = gridy_max - gridy_min

    marknum = model.markers.parameters.distribution.marknum.value
    for mm1 in 1:marknum
        y_marker = marker_y[mm1]
        if y_marker > gridy_max
            marker_y[mm1] = calculate_y(marker_y[mm1], ysize)
            marker_temperature[mm1] = calculate_temperature(
                marker_temperature[mm1], t_top, t_bottom)
        end
    end
    return nothing
end

function calculate_y(y_marker::Float64, ysize::Float64)::Float64
    return y_marker - ysize
end

function calculate_temperature(
    temperature::Float64,
    temperature_top::Float64,
    temperature_bottom::Float64
)::Float64
    return temperature + temperature_top - temperature_bottom
end

end # module Channel 