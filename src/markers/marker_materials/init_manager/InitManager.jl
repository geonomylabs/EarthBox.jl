module InitManager

include("utils/__init__.jl")
include("types/LithosphericExtensionWeakFault.jl")
include("types/LithosphericExtensionWeakSeed.jl")
include("types/LithosphericExtensionNoSeed.jl")
include("types/LithosphericExtensionAndPlume.jl")
include("types/LithosphericExtension.jl")
include("types/FlexureTriangularHole.jl")
include("types/PlasticityBenchmarkKaus10.jl")
include("types/ElasticSlab.jl")
include("types/FractureZone.jl")
include("types/RayleighTaylorInstability.jl")
include("types/SandboxExtension.jl")
include("types/SandboxShortening.jl")
include("types/SeafloorSpreading.jl")
include("types/SimplePlume.jl")
include("types/SimpleSediment.jl")
include("types/StickyOnMantle.jl")
include("types/Uniform.jl")
include("types/ViscousBlock.jl")

import EarthBox.ModelDataContainer: ModelData
import ..Options: option_names

function initialize!(model::ModelData, ::Val{option_names.LithosphericExtensionWeakFault})
    LithosphericExtensionWeakFault.initialize!(model)
end

function initialize!(model::ModelData, ::Val{option_names.LithosphericExtensionWeakSeed})
    LithosphericExtensionWeakSeed.initialize!(model)
end

function initialize!(model::ModelData, ::Val{option_names.LithosphericExtensionNoSeed})
    LithosphericExtensionNoSeed.initialize!(model)
end

function initialize!(model::ModelData, ::Val{option_names.LithosphericExtensionMantlePlume})
    LithosphericExtensionAndPlume.initialize!(model)
end

function initialize!(model::ModelData, ::Val{option_names.LithosphericExtensionLateralStrongZones})
    LithosphericExtension.initialize!(model)
end

function initialize!(model::ModelData, ::Val{option_names.FlexureTriangularHole})
    FlexureTriangularHole.initialize!(model)
end

function initialize!(model::ModelData, ::Val{option_names.PlasticityBenchmarkWeakNotchKaus10})
    PlasticityBenchmarkKaus10.initialize!(model)
end

function initialize!(model::ModelData, ::Val{option_names.ElasticSlab})
    ElasticSlab.initialize!(model)
end

function initialize!(model::ModelData, ::Val{option_names.FractureZone})
    FractureZone.initialize!(model)
end

function initialize!(model::ModelData, ::Val{option_names.RayleighTaylorInstabilityBenchmark})
    RayleighTaylorInstability.initialize!(model)
end

function initialize!(model::ModelData, ::Val{option_names.SandboxExtension})
    SandboxExtension.initialize!(model)
end

function initialize!(model::ModelData, ::Val{option_names.SandboxShortening})
    SandboxShortening.initialize!(model)
end

function initialize!(model::ModelData, ::Val{option_names.SeafloorSpreading})
    SeafloorSpreading.initialize!(model)
end

function initialize!(model::ModelData, ::Val{option_names.SimplePlume})
    SimplePlume.initialize!(model)
end

function initialize!(model::ModelData, ::Val{option_names.SimpleSediment})
    SimpleSediment.initialize!(model)
end

function initialize!(model::ModelData, ::Val{option_names.StickyOnMantle})
    StickyOnMantle.initialize!(model)
end

function initialize!(model::ModelData, ::Val{option_names.UniformComposition})
    Uniform.initialize!(model)
end

end # module