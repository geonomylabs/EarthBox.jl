module MultipleTimeStepIncrease

import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: load_parameters!

struct ValidInputNames
    ntime_increase_1::Symbol
    ntime_increase_2::Symbol
    ntime_increase_3::Symbol
    ntime_increase_4::Symbol
    time_increase_factor::Symbol
    iuse_multiple_timestep_increase::Symbol
end

function initialize!(model::ModelData, kwargs...)
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
end

end