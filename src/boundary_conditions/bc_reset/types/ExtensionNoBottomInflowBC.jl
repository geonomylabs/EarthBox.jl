module ExtensionNoBottomInflowBC

import EarthBox.ModelDataContainer: ModelData
import ...SetBoundaryConditions
import ...BoundaryVelocity: ExtensionFreeSlipAlongBottom

"""
    set_bcs!(model::ModelData)

Boundary condition coefficients for extension with a notch.

Simple extension notch experiment with prescribed side and top velocity. Free 
slip on bottom boundary. Velocity along top boundary ensures volume conservation.
"""
function set_bcs!(model::ModelData)
    vx = model.bcs.parameters.velocity.full_velocity_extension.value/2.0

    # Calculate velocity along top boundary to ensure volume conservation
    ExtensionFreeSlipAlongBottom.calc_velocity!(model)
    vyu = model.bcs.parameters.velocity.vyu.value

    SetBoundaryConditions.Pressure.upper_left_corner!(model)
    SetBoundaryConditions.TopBoundary.inflow_outflow!(model, vyu)
    SetBoundaryConditions.BottomBoundary.free_slip!(model)
    SetBoundaryConditions.LeftBoundary.inflow_outflow!(model, -vx)
    SetBoundaryConditions.RightBoundary.inflow_outflow!(model, vx)
    SetBoundaryConditions.TopBoundary.temperature!(model)
    SetBoundaryConditions.BottomBoundary.temperature!(model)
    SetBoundaryConditions.LeftBoundary.zero_heat_flux!(model)
    SetBoundaryConditions.RightBoundary.zero_heat_flux!(model)
end

end # module ExtensionNoBottomInflowBC 