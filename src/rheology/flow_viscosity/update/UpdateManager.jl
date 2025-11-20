module UpdateManager

include("utils/FlowUtils.jl")
include("types/blankenbach89/Blankenback89.jl")
include("types/couette_flow_benchmark/CouetteFlowBenchmark.jl")
include("types/isoviscous/Isoviscous.jl")
include("types/x_dependent/Xdependent.jl")
include("types/composite_viscosity/CompositeViscosity.jl")

import Profile
import EarthBox.ModelDataContainer: ModelData
import ..Options: option_names
import .Isoviscous
import .CompositeViscosity
import .Blankenback89
import .CouetteFlowBenchmark
import .Xdependent

function update_marker_flow_viscosity!(
    model::ModelData, 
    inside_flags::Vector{Int8},
    ::Val{option_names.Isoviscous}
)::Nothing
    Isoviscous.update_marker_flow_viscosity!(model, inside_flags)
    return nothing
end

function update_marker_flow_viscosity!(
    model::ModelData, 
    inside_flags::Vector{Int8},
    ::Val{option_names.Composite}
)::Nothing
    #Profile.@profile 
    CompositeViscosity.update_marker_flow_viscosity!(model, inside_flags)
    #Profile.print(format=:flat, C=true, combine=true, sortedby=:count)
    return nothing
end

function update_marker_flow_viscosity!(
    model::ModelData, 
    inside_flags::Vector{Int8},
    ::Val{option_names.TemperatureDependentConvectionBenchmark}
)::Nothing
    Blankenback89.update_marker_flow_viscosity!(model, inside_flags)
    return nothing
end

function update_marker_flow_viscosity!(
    model::ModelData, 
    inside_flags::Vector{Int8},
    ::Val{option_names.TemperatureDependentCouetteBenchmark}
)::Nothing
    CouetteFlowBenchmark.update_marker_flow_viscosity!(model, inside_flags)
    return nothing
end

function update_marker_flow_viscosity!(
    model::ModelData, 
    inside_flags::Vector{Int8},
    ::Val{option_names.XCoordinateDependent}
)::Nothing
    Xdependent.update_marker_flow_viscosity!(model, inside_flags)
    return nothing
end

end # module

