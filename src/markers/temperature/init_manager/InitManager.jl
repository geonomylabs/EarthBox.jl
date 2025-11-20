module InitManager

include("utils/TempIniStructs.jl")
include("utils/PlumeUpdate.jl")
include("seafloor_spreading/HalfSpaceCooling.jl")
include("steady_state/GeothermStructs.jl")
include("steady_state/Geotherm.jl")
include("types/TemperatureWave.jl")
include("types/Uniform.jl")
include("types/BoxConvection.jl")
include("types/Linear.jl")
include("types/AnalyticalThreeLayer.jl")
include("types/FourLinearSegments.jl")
include("types/FractureZone.jl")
include("types/HotBox.jl")
include("types/HalfSpaceSpreading.jl")
include("types/HalfSpaceDoubleSpreading.jl")

import EarthBox.ModelDataContainer: ModelData
import ..Options: option_names

function initialize!(model::ModelData, ::Val{option_names.TemperatureWave})
    TemperatureWave.initialize!(model)
end

function initialize!(model::ModelData, ::Val{option_names.Uniform})
    Uniform.initialize!(model)
end

function initialize!(model::ModelData, ::Val{option_names.BoxConvection})
    BoxConvection.initialize!(model)
end

function initialize!(model::ModelData, ::Val{option_names.Linear})
    Linear.initialize!(model)
end

function initialize!(model::ModelData, ::Val{option_names.AnalyticalThreeLayer})
    AnalyticalThreeLayer.initialize!(model)
end

function initialize!(model::ModelData, ::Val{option_names.FourLinearSegments})
    FourLinearSegments.initialize!(model)
end

function initialize!(model::ModelData, ::Val{option_names.FractureZone})
    FractureZone.initialize!(model)
end

function initialize!(model::ModelData, ::Val{option_names.HotBox})
    HotBox.initialize!(model)
end

function initialize!(model::ModelData, ::Val{option_names.HalfSpaceSpreading})
    HalfSpaceSpreading.initialize!(model)
end

function initialize!(model::ModelData, ::Val{option_names.HalfSpaceDoubleSpreading})
    HalfSpaceDoubleSpreading.initialize!(model)
end

end








