module LithosphericExtensionFixedBCAsymmetric

import EarthBox.ModelDataContainer: ModelData
import ...SetBoundaryConditions
import ...BoundaryVelocity: ExtensionInflowAlongBottom

function set_bcs!(model::ModelData)
    vx = model.bcs.parameters.velocity.full_velocity_extension.value
    ExtensionInflowAlongBottom.calc_velocity!(model, use_asymmetric_inflow_outflow=true)
    vyu = model.bcs.parameters.velocity.vyu.value
    vyl = model.bcs.parameters.velocity.vyl.value
    SetBoundaryConditions.Pressure.upper_left_corner!(model)
    SetBoundaryConditions.TopBoundary.inflow_outflow!(model, vyu)
    SetBoundaryConditions.BottomBoundary.inflow_outflow!(model, vyl)
    SetBoundaryConditions.LeftBoundary.inflow_outflow!(model, 0.0)
    SetBoundaryConditions.RightBoundary.inflow_outflow!(model, vx)
    SetBoundaryConditions.TopBoundary.temperature!(model)
    SetBoundaryConditions.BottomBoundary.temperature!(model)
    SetBoundaryConditions.LeftBoundary.zero_heat_flux!(model)
    SetBoundaryConditions.RightBoundary.zero_heat_flux!(model)
end

end # module LithosphericExtensionFixedBCAsymmetric 