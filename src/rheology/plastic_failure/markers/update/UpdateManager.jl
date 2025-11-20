module UpdateManager

include("utils/Limits.jl")
include("utils/YieldCheck.jl")
include("types/ElasticForecast.jl")
include("types/ViscoelasticForecast.jl")

import EarthBox.ModelDataContainer: ModelData
import ..Options: option_names

function update_plastic_yielding!(
    model::ModelData,
    inside_flags::Vector{Int8},
    no_yielding_in_mobile_wall::Bool,
    ::Val{option_names.PureElasticStressForecast}
)::Nothing
    ElasticForecast.update_plastic_yielding!(model, inside_flags, no_yielding_in_mobile_wall)
    return nothing
end

function update_plastic_yielding!(
    model::ModelData,
    inside_flags::Vector{Int8},
    no_yielding_in_mobile_wall::Bool,
    ::Val{option_names.ViscoelasticStressForecast}
)::Nothing
    ViscoelasticForecast.update_plastic_yielding!(model, inside_flags, no_yielding_in_mobile_wall)
    return nothing
end

end # module
