module LithosphericExtensionInflowAndOutflowAlongSidesBC

import EarthBox.ModelDataContainer: ModelData
import ...SetBoundaryConditions
import ...BoundaryVelocity: ExtensionInflowAndOutflowAlongSides

function set_bcs!(model::ModelData)
    velocity_extension = abs(model.bcs.parameters.velocity.full_velocity_extension.value)/2.0
    ExtensionInflowAndOutflowAlongSides.calc_velocity!(model)
    (
        velocity_y_upper, velocity_inflow_left, velocity_inflow_right
    ) = get_conservative_velocities(model)
    SetBoundaryConditions.Pressure.upper_left_corner!(model)
    SetBoundaryConditions.TopBoundary.inflow_outflow!(model, velocity_y_upper)
    SetBoundaryConditions.BottomBoundary.free_slip!(model)
    SetBoundaryConditions.LeftBoundary.inflow_and_outflow_along_side!(
        model, -velocity_extension, velocity_inflow_left)
    SetBoundaryConditions.RightBoundary.inflow_and_outflow_along_side!(
        model, velocity_extension, velocity_inflow_right)
    SetBoundaryConditions.TopBoundary.temperature!(model)
    SetBoundaryConditions.BottomBoundary.temperature!(model)
    SetBoundaryConditions.LeftBoundary.zero_heat_flux!(model)
    SetBoundaryConditions.RightBoundary.zero_heat_flux!(model)
end

function get_conservative_velocities(model::ModelData)::Tuple{Float64, Float64, Float64}
    velocity = model.bcs.parameters.velocity
    velocity_y_upper = velocity.vyu.value
    velocity_inflow_left = velocity.velocity_inflow_left.value
    velocity_inflow_right = velocity.velocity_inflow_right.value
    return velocity_y_upper, velocity_inflow_left, velocity_inflow_right
end

end # module