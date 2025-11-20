module BCResetManager

include("types/BoxConvectionBC.jl")
include("types/ElasticSlabBC.jl")
include("types/PureShearDeformationBC.jl")
include("types/RayleighTaylorInstabilityBC.jl")
include("types/SandboxBC.jl")
include("types/SolidBodyRotationBC.jl")
include("types/ExtensionNoBottomInflowBC.jl")
include("types/ViscousBlockBC.jl")
include("types/VerticalChannelFlowIsoviscousBC.jl")
include("types/VerticalChannelFlowNonNewtonianBC.jl")
include("types/VerticalCouetteFlowBC.jl")
include("types/VerticalChannelFlowShearHeatingBC.jl")
include("types/LithosphericContractionFixedBCAsymmetric.jl")
include("types/LithosphericContractionFixedBC.jl")
include("types/LithosphericExtensionFixedBC.jl")
include("types/LithosphericExtensionFixedBCAsymmetric.jl")
include("types/LithosphericExtensionDepthDependentBC.jl")
include("types/LithosphericExtensionInflowAndOutflowAlongSidesBC.jl")
include("types/LithosphericExtensionMovingBC.jl")

import EarthBox.ModelDataContainer: ModelData
import ..Options: option_names

function reset_bcs!(model::ModelData, ::Val{option_names.LithosphericExtensionMovingBoundaries})::Nothing
    LithosphericExtensionMovingBC.set_bcs!(model)
    return nothing
end

function reset_bcs!(model::ModelData, ::Val{option_names.LithosphericExtensionFixedBoundaries})::Nothing
    LithosphericExtensionFixedBC.set_bcs!(model)
    return nothing
end

function reset_bcs!(model::ModelData, ::Val{option_names.LithosphericExtensionFixedBoundariesAsymmetric})::Nothing
    LithosphericExtensionFixedBCAsymmetric.set_bcs!(model)
    return nothing
end

function reset_bcs!(model::ModelData, ::Val{option_names.LithosphericExtensionDepthDependent})::Nothing
    LithosphericExtensionDepthDependentBC.set_bcs!(model)
    return nothing
end

function reset_bcs!(model::ModelData, ::Val{option_names.LithosphericExtensionInflowAndOutflowAlongSides})::Nothing
    LithosphericExtensionInflowAndOutflowAlongSidesBC.set_bcs!(model)
    return nothing
end

function reset_bcs!(model::ModelData, ::Val{option_names.ExtensionNoBottomInflow})::Nothing
    ExtensionNoBottomInflowBC.set_bcs!(model)
    return nothing
end

function reset_bcs!(model::ModelData, ::Val{option_names.PureShearDeformation})::Nothing
    PureShearDeformationBC.set_bcs!(model)
    return nothing
end

function reset_bcs!(model::ModelData, ::Val{option_names.VerticalCouetteFlow})::Nothing
    VerticalCouetteFlowBC.set_bcs!(model)
    return nothing
end

function reset_bcs!(model::ModelData, ::Val{option_names.SolidBodyRotation})::Nothing
    SolidBodyRotationBC.set_bcs!(model)
    return nothing
end

function reset_bcs!(model::ModelData, ::Val{option_names.RayleighTaylorInstability})::Nothing
    RayleighTaylorInstabilityBC.set_bcs!(model)
    return nothing
end

function reset_bcs!(model::ModelData, ::Val{option_names.BoxConvection})::Nothing
    BoxConvectionBC.set_bcs!(model)
    return nothing
end

function reset_bcs!(model::ModelData, ::Val{option_names.ElasticSlab})::Nothing
    ElasticSlabBC.set_bcs!(model)
    return nothing
end

function reset_bcs!(model::ModelData, ::Val{option_names.ViscousBlock})::Nothing
    ViscousBlockBC.set_bcs!(model)
    return nothing
end

function reset_bcs!(model::ModelData, ::Val{option_names.SandboxShortening})::Nothing
    SandboxBC.set_bcs!(model)
    return nothing
end

function reset_bcs!(model::ModelData, ::Val{option_names.SandboxExtension})::Nothing
    SandboxBC.set_bcs!(model)
    return nothing
end

function reset_bcs!(model::ModelData, ::Val{option_names.VerticalChannelFlowIsoviscous})::Nothing
    VerticalChannelFlowIsoviscousBC.set_bcs!(model)
    return nothing
end

function reset_bcs!(model::ModelData, ::Val{option_names.VerticalChannelFlowNonNewtonian})::Nothing
    VerticalChannelFlowNonNewtonianBC.set_bcs!(model)
    return nothing
end

function reset_bcs!(model::ModelData, ::Val{option_names.VerticalChannelFlowShearHeating})::Nothing
    VerticalChannelFlowShearHeatingBC.set_bcs!(model)
    return nothing
end

function reset_bcs!(model::ModelData, ::Val{option_names.FractureZoneContractionFixedBoundaries})::Nothing
    LithosphericContractionFixedBC.set_bcs!(model)
    return nothing
end

function reset_bcs!(model::ModelData, ::Val{option_names.FractureZoneContractionFixedBoundariesAsymmetric})::Nothing
    LithosphericContractionFixedBCAsymmetric.set_bcs!(model)
    return nothing
end

function reset_bcs!(model::ModelData, ::Val{option_names.ViscoelasticContraction})::Nothing
    LithosphericContractionFixedBC.set_bcs!(model)
    return nothing
end

function reset_bcs!(model::ModelData, ::Val{option_names.ViscoelasticContractionAsymmetric})::Nothing
    LithosphericContractionFixedBCAsymmetric.set_bcs!(model)
    return nothing
end

end # module

