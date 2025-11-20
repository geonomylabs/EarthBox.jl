module SingleTimeStepIncrease

import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: load_parameters!

struct ValidInputNames
    iuse_single_timestep_increase::Symbol
    ntime_increase::Symbol
    time_increase_factor::Symbol
end

function initialize!(model::ModelData; kwargs...)
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
end

end # module 