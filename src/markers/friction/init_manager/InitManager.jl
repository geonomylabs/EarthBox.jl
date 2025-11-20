module InitManager

include("types/Randomized.jl")
include("types/Regular.jl")

import EarthBox.ModelDataContainer: ModelData
import ..Options: option_names
import .Randomized
import .Regular

function initialize!(
    model::ModelData, 
    ::Val{option_names.Randomized}
)
    Randomized.initialize!(model)
end

function initialize!(
    model::ModelData, 
    ::Val{option_names.Regular}
)
    Regular.initialize!(model)
end

end # module
