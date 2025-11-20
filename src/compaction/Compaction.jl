module Compaction

include("core/CompactionTools.jl")
include("core/MarkerCompaction.jl")
include("core/UpdateProperties.jl")
include("core/CompactionCorrection.jl")
include("core/ApplyCompaction.jl")

import EarthBox.ModelDataContainer: ModelData

function get_boolean_options_for_compaction(model::ModelData)
    downhill_diffusion = model.topography.parameters.downhill_diffusion
    use_compaction_correction = Bool(downhill_diffusion.iuse_compaction_correction.value)
    return use_compaction_correction
end

end