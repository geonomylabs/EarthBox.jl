module LoopInputStruct

import EarthBox.ConfigurationManager.SolverConfig: SolverConfigState

struct LoopInput
    run_paths::Dict{String, String}
    solver_config::SolverConfigState
    no_yielding_in_mobile_wall::Bool
    no_yielding_in_plate_extension::Bool
end

end

