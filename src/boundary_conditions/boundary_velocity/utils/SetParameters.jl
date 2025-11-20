module SetParameters

import EarthBox.ModelDataContainer: ModelData

"""
    set_velocity_upper_and_lower(model::ModelData, velocity_y_upper::Float64, velocity_y_lower::Float64)

Set velocity at top and bottom boundaries of extension model.

# Arguments
- `model::ModelData`: Model data container
- `velocity_y_upper::Float64`:
  - Conservative vertical velocity (m/s) along upper boundary of extension model
- `velocity_y_lower::Float64`:
  - Conservative vertical velocity (m/s) along lower boundary of extension model

# Updated parameter from group `bcs.parameters.velocity`
- `vyu.value::Float64`:
  - Conservative vertical velocity (m/s) along upper boundary of extension model
- `vyl.value::Float64`:
  - Conservative vertical velocity (m/s) along lower boundary of extension model
"""
function set_velocity_upper_and_lower(
    model::ModelData,
    velocity_y_upper::Float64,
    velocity_y_lower::Float64
)
    model.bcs.parameters.velocity.vyu.value = velocity_y_upper
    model.bcs.parameters.velocity.vyl.value = velocity_y_lower
end

"""
    set_velocity_upper(model::ModelData, velocity_y_upper::Float64)

Set velocity at top and bottom boundaries of extension model.

# Arguments
- `model::ModelData`: Model data container
- `velocity_y_upper::Float64`:
  - Conservative vertical velocity (m/s) along upper boundary of extension model

# Updated parameter from group `bcs.parameters.velocity`
- `vyu.value::Float64`:
  - Conservative vertical velocity (m/s) along upper boundary of extension model
"""
function set_velocity_upper(
    model::ModelData,
    velocity_y_upper::Float64
)
    model.bcs.parameters.velocity.vyu.value = velocity_y_upper
end

"""
    set_inflow_velocity(
        model::ModelData, velocity_inflow_left::Float64, velocity_inflow_right::Float64)

Set velocity at top and bottom boundaries of extension model.

# Arguments
- `model::ModelData`: Model data container
- `velocity_inflow_left::Float64`: Inflow velocity along left boundary (m/s)
- `velocity_inflow_right::Float64`: Inflow velocity along right boundary (m/s)

# Updated parameter from group `bcs.parameters.velocity`
- `velocity_inflow_left.value::Float64`: Inflow velocity along left boundary (m/s)
- `velocity_inflow_right.value::Float64`: Inflow velocity along right boundary (m/s)
"""
function set_inflow_velocity(
    model::ModelData,
    velocity_inflow_left::Float64,
    velocity_inflow_right::Float64
)
    model.bcs.parameters.velocity.velocity_inflow_left.value = velocity_inflow_left
    model.bcs.parameters.velocity.velocity_inflow_right.value = velocity_inflow_right
end

end # module SetParameters 