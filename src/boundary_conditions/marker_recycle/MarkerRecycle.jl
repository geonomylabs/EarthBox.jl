module MarkerRecycle

include("utils/CheckDicts.jl")
include("utils/RandomMarkerArray.jl")
include("utils/PressureReset.jl")
include("utils/Sticky.jl")
include("utils/Coordinates.jl")
include("utils/CountRecycle.jl")
include("utils/ContractionRecycleLocations.jl")
include("utils/ResetMarkers.jl")
include("utils/FractureZoneProperties.jl")
include("utils/InitializeMarkerRecycling.jl")
include("types/Channel.jl")
include("types/Extension.jl")
include("types/ExtensionNoBottomInflow.jl")
include("types/ExtensionInflowAndOutflowAlongSides.jl")
include("types/ContractionFractureZone.jl")
include("types/ContractionViscoelastic.jl")

import EarthBox.ModelDataContainer: ModelData
import ..Options: option_names

function recycle_markers!(
    model::ModelData, 
    inside_flags::Vector{Int8}, 
    ::Val{option_names.LithosphericExtensionMovingBoundaries}
)::Nothing
    return nothing
end

function recycle_markers!(
    model::ModelData, 
    inside_flags::Vector{Int8}, 
    ::Val{option_names.LithosphericExtensionFixedBoundaries}
)::Nothing
    Extension.recycle_markers!(model, inside_flags)
    return nothing
end

function recycle_markers!(
    model::ModelData, 
    inside_flags::Vector{Int8}, 
    ::Val{option_names.LithosphericExtensionFixedBoundariesAsymmetric}
)::Nothing
    Extension.recycle_markers!(model, inside_flags)
    return nothing
end

function recycle_markers!(
    model::ModelData, 
    inside_flags::Vector{Int8}, 
    ::Val{option_names.LithosphericExtensionDepthDependent}
)::Nothing
    Extension.recycle_markers!(model, inside_flags)
    return nothing
end

function recycle_markers!(
    model::ModelData, 
    inside_flags::Vector{Int8}, 
    ::Val{option_names.LithosphericExtensionInflowAndOutflowAlongSides}
)::Nothing
    ExtensionInflowAndOutflowAlongSides.recycle_markers!(model, inside_flags)
    return nothing
end

function recycle_markers!(
    model::ModelData, 
    inside_flags::Vector{Int8}, 
    ::Val{option_names.ExtensionNoBottomInflow}
)::Nothing
    ExtensionNoBottomInflow.recycle_markers!(model, inside_flags)
    return nothing
end

function recycle_markers!(
    model::ModelData, 
    inside_flags::Vector{Int8}, 
    ::Val{option_names.PureShearDeformation}
)::Nothing
    return nothing
end

function recycle_markers!(
    model::ModelData, 
    inside_flags::Vector{Int8}, 
    ::Val{option_names.VerticalCouetteFlow}
)::Nothing
    Channel.recycle_markers!(model)
    return nothing
end

function recycle_markers!(
    model::ModelData, 
    inside_flags::Vector{Int8}, 
    ::Val{option_names.SolidBodyRotation}
)::Nothing
    return nothing
end

function recycle_markers!(
    model::ModelData,
    inside_flags::Vector{Int8}, 
    ::Val{option_names.RayleighTaylorInstability}
)::Nothing
    return nothing
end

function recycle_markers!(
    model::ModelData, 
    inside_flags::Vector{Int8}, 
    ::Val{option_names.BoxConvection}
)::Nothing
    return nothing
end

function recycle_markers!(
    model::ModelData, 
    inside_flags::Vector{Int8}, 
    ::Val{option_names.ElasticSlab}
)::Nothing
    return nothing
end

function recycle_markers!(
    model::ModelData,
    inside_flags::Vector{Int8}, 
    ::Val{option_names.ViscousBlock}
)::Nothing
    return nothing
end

function recycle_markers!(
    model::ModelData, 
    inside_flags::Vector{Int8}, 
    ::Val{option_names.SandboxShortening}
)::Nothing
    return nothing
end

function recycle_markers!(
    model::ModelData, 
    inside_flags::Vector{Int8}, 
    ::Val{option_names.SandboxExtension}
)::Nothing
    return nothing
end

function recycle_markers!(
    model::ModelData,
    inside_flags::Vector{Int8}, 
    ::Val{option_names.VerticalChannelFlowIsoviscous}
)::Nothing
    Channel.recycle_markers!(model)
    return nothing
end

function recycle_markers!(
    model::ModelData, 
    inside_flags::Vector{Int8}, 
    ::Val{option_names.VerticalChannelFlowNonNewtonian}
)::Nothing
    Channel.recycle_markers!(model)
    return nothing
end

function recycle_markers!(
    model::ModelData, 
    inside_flags::Vector{Int8}, 
    ::Val{option_names.VerticalChannelFlowShearHeating}
)::Nothing
    Channel.recycle_markers!(model)
    return nothing
end

function recycle_markers!(
    model::ModelData, 
    inside_flags::Vector{Int8}, 
    ::Val{option_names.FractureZoneContractionFixedBoundaries}
)::Nothing
    ContractionFractureZone.recycle_markers_symmetric!(model, inside_flags)
    return nothing
end

function recycle_markers!(
    model::ModelData, 
    inside_flags::Vector{Int8}, 
    ::Val{option_names.FractureZoneContractionFixedBoundariesAsymmetric}
)::Nothing
    ContractionFractureZone.recycle_markers_asymmetric!(model, inside_flags)
    return nothing
end

function recycle_markers!(
    model::ModelData, 
    inside_flags::Vector{Int8}, 
    ::Val{option_names.ViscoelasticContraction}
)::Nothing
    ContractionViscoelastic.recycle_markers_symmetric!(model, inside_flags)
    return nothing
end

function recycle_markers!(
    model::ModelData, 
    inside_flags::Vector{Int8}, 
    ::Val{option_names.ViscoelasticContractionAsymmetric}
)::Nothing
    ContractionViscoelastic.recycle_markers_asymmetric!(model, inside_flags)
    return nothing
end

end

