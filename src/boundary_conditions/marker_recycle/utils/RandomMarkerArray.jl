module RandomMarkerArray

import EarthBox.ModelDataContainer: ModelData

function get_random_marker_array(model::ModelData)::Vector{Float64}
    marker_params = model.markers.parameters
    mxnum = marker_params.distribution.mxnum.value
    mynum = marker_params.distribution.mynum.value
    return rand(mxnum * mynum)
end

end # module 