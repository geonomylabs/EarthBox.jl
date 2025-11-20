module MarkersRandomized

import EarthBox.ModelDataContainer: ModelData

function initialize!(model::ModelData, marker_random::Vector{Float64})
    mxnum = model.markers.parameters.distribution.mxnum.value
    mynum = model.markers.parameters.distribution.mynum.value

    mxstep = model.markers.parameters.distribution.mxstep.value
    mystep = model.markers.parameters.distribution.mystep.value
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array

    imarker = 1
    for ixmarker in 1:mxnum
        for iymarker in 1:mynum
            random_value = marker_random[imarker]
            marker_x[imarker] = calc_coordinate_random(
                ixmarker, mxstep, random_value)
            marker_y[imarker] = calc_coordinate_random(
                iymarker, mystep, random_value)
            imarker += 1
        end
    end
end

function calc_coordinate_random(
    imarker::Int,
    mstep::Float64,
    random_value::Float64
)::Float64
    coor = (
        Float64(imarker)*mstep
        - mstep / 2.0
        + (random_value - 0.5) * mstep
    )
    return coor
end

end # module 