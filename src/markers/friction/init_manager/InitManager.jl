module InitManager

include("types/Randomized.jl")
include("types/Regular.jl")
include("types/CentralWeakening.jl")

import EarthBox.ModelDataContainer: ModelData
import ..Options: option_names
import .Randomized
import .Regular
import .CentralWeakening

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

function initialize!(
    model::ModelData, 
    ::Val{option_names.CentralWeakening}
)
    CentralWeakening.initialize!(model)
end

end # module
