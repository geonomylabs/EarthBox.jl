module Internal

import EarthBox.ModelDataContainer: ModelData
import ..InternalBC

function deactivate!(model::ModelData)::Nothing
    InternalBC.turn_off_internal_velocity!(model)
    return nothing
end

function mobile_wall_sandbox!(model::ModelData)::Nothing
    InternalBC.set_internal_velocity_for_sandbox_with_mobile_wall!(model)
    return nothing
end

end # module Internal 