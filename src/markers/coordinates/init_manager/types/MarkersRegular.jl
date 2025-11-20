module MarkersRegular

import EarthBox.ModelDataContainer: ModelData

function initialize!(model::ModelData)
    mxnum = model.markers.parameters.distribution.mxnum.value
    mynum = model.markers.parameters.distribution.mynum.value

    mxstep = model.markers.parameters.distribution.mxstep.value
    mystep = model.markers.parameters.distribution.mystep.value
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array

    imarker = 1  # Julia uses 1-based indexing
    for ixmarker in 1:mxnum
        for iymarker in 1:mynum
            marker_x[imarker] = calc_coordinate_regular(ixmarker, mxstep)
            marker_y[imarker] = calc_coordinate_regular(iymarker, mystep)
            imarker += 1
        end
    end
end

function calc_coordinate_regular(imarker::Int, mstep::Float64)::Float64
    coor = Float64(imarker) * mstep - mstep / 2.0
    return coor
end

end # module 