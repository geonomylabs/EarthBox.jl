module LithosphericExtensionMovingBC

import EarthBox.ModelDataContainer: ModelData
import ...SetBoundaryConditions

function set_bcs!(model::ModelData)
    vx = model.bcs.parameters.velocity.full_velocity_extension.value
    vyu = model.bcs.parameters.velocity.vyu.value
    vyl = model.bcs.parameters.velocity.vyl.value
    SetBoundaryConditions.Pressure.upper_left_corner!(model)
    SetBoundaryConditions.TopBoundary.inflow_outflow!(model, vyu)
    SetBoundaryConditions.BottomBoundary.inflow_outflow!(model, vyl)
    SetBoundaryConditions.LeftBoundary.inflow_outflow!(model, -vx/2.0)
    SetBoundaryConditions.RightBoundary.inflow_outflow!(model, vx/2.0)
    SetBoundaryConditions.TopBoundary.temperature!(model)
    SetBoundaryConditions.BottomBoundary.temperature!(model)
    SetBoundaryConditions.LeftBoundary.zero_heat_flux!(model)
    SetBoundaryConditions.RightBoundary.zero_heat_flux!(model)
end

end # module