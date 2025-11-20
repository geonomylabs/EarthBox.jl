module Pressure

import EarthBox.ModelDataContainer: ModelData
import ..PressureBC

function upper_left_corner!(model::ModelData)::Nothing
    PressureBC.set_pressure_bc_mode!(model, 0)
    return nothing
end

function top_and_bottom!(model::ModelData)::Nothing
    PressureBC.set_pressure_bc_mode!(model, 1)
    return nothing
end

function left_and_right!(model::ModelData)::Nothing
    PressureBC.set_pressure_bc_mode!(model, 2)
    return nothing
end

end # module Pressure 