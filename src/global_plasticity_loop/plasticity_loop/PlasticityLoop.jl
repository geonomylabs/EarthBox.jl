module PlasticityLoop

include("plasticity_loop_markers/PlasticityLoopMarkers.jl")
include("plasticity_loop_nodes/PlasticityLoopNodes.jl")

import EarthBox.ModelDataContainer: ModelData
import EarthBox.LoopInputStruct: LoopInput
import ..Options: option_names

function run_picard_loop(
    model::ModelData, 
    loop_input::LoopInput, 
    inside_flags::Vector{Int8},
    ::Val{option_names.MarkerPlasticityLoop}
)::Nothing
    PlasticityLoopMarkers.picard_loop(model, loop_input, inside_flags)
    return nothing
end

function run_picard_loop(
    model::ModelData, 
    loop_input::LoopInput,
    inside_flags::Vector{Int8},
    ::Val{option_names.NodalPlasticityLoop}
)::Nothing
    PlasticityLoopNodes.picard_loop(model, loop_input, inside_flags)
    return nothing
end

end

