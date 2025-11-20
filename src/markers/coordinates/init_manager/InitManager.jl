module InitManager

include("types/MarkersRandomized.jl")
include("types/MarkersRegular.jl")

import Random: rand
import EarthBox.ModelDataContainer: ModelData
import ..Options: option_names
import .MarkersRandomized
import .MarkersRegular

function initialize!(
    model::ModelData,
    ::Val{option_names.Randomized}
)
    mxnum = model.markers.parameters.distribution.mxnum.value
    mynum = model.markers.parameters.distribution.mynum.value
    marker_random = rand(mxnum*mynum)
    MarkersRandomized.initialize!(model, marker_random)   
end

function initialize!(
    model::ModelData,
    ::Val{option_names.Regular}
)
    MarkersRegular.initialize!(model)
end

end # module
