module Markers

include("utils/RayleighTaylor.jl")
include("utils/Layers.jl")
include("coordinates/MarkerCoordinates.jl")
include("marker_materials/MarkerMaterials.jl")
include("friction/MarkerFriction.jl")
include("temperature/MarkerTemperature.jl")
include("cohesion/MarkerCohesion.jl")
include("dilatancy/MarkerDilatancy.jl")
include("preexponential/MarkerPreexponential.jl")

import EarthBox.ModelDataContainer: ModelData
import .MarkerCoordinates
import .MarkerMaterials
import .MarkerMaterials: MarkerStressLimits
import .MarkerMaterials: MarkerBoundaryFriction
import .MarkerMaterials: MarkerViscousStrainSoftening
import .MarkerFriction
import .MarkerTemperature
import .MarkerCohesion
import .MarkerDilatancy
import .MarkerPreexponential

function initialize!(
    model::ModelData, 
    paths::Dict{String, String}
)::Nothing
    MarkerCoordinates.initialize!(model)
    MarkerMaterials.initialize!(model, paths=paths)
    MarkerStressLimits.initialize!(model)
    MarkerBoundaryFriction.initialize!(model)
    MarkerFriction.initialize!(model)
    MarkerViscousStrainSoftening.initialize!(model)
    MarkerTemperature.initialize!(model)
    MarkerCohesion.initialize!(model)
    MarkerDilatancy.initialize!(model)
    MarkerPreexponential.initialize!(model)
    nothing
end

end # module Markers 