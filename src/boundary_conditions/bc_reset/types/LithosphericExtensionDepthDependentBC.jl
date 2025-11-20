module LithosphericExtensionDepthDependentBC

import EarthBox.ModelDataContainer: ModelData
import ...SetBoundaryConditions
import ...BoundaryVelocity: ExtensionDepthDependentWithInflow

"""
    set_bcs!(model::ModelData)

Boundary condition coefficients for depth-dependent lithospheric extension with 
fixed boundaries.

Extension models with fixed boundaries. Side boundaries have free slip with
depth-dependent orthogonal velocity prescribed. Along top and bottom
boundaries free slip is applied with prescribed orthogonal velocity
calculated to ensure volume conservation. Particles are recycled (no
regridding) with special treatment in sticky air region to ensure constant
sticky air volume. Temperature boundary conditions on top and bottom, zero flux
on sides.
"""
function set_bcs!(model::ModelData)
    vx = model.bcs.parameters.velocity.full_velocity_extension.value
    y_linear_left, y_linear_right = ExtensionDepthDependentWithInflow.calc_velocity!(model)
    vyu = model.bcs.parameters.velocity.vyu.value
    vyl = model.bcs.parameters.velocity.vyl.value
    SetBoundaryConditions.Pressure.upper_left_corner!(model)
    SetBoundaryConditions.TopBoundary.inflow_outflow!(model, vyu)
    SetBoundaryConditions.BottomBoundary.inflow_outflow!(model, vyl)
    SetBoundaryConditions.LeftBoundary.inflow_outflow_depth_dependent!(model, y_linear_left, vx/2.0)
    SetBoundaryConditions.RightBoundary.inflow_outflow_depth_dependent!(model, y_linear_right, vx/2.0)
    SetBoundaryConditions.TopBoundary.temperature!(model)
    SetBoundaryConditions.BottomBoundary.temperature!(model)
    SetBoundaryConditions.LeftBoundary.zero_heat_flux!(model)
    SetBoundaryConditions.RightBoundary.zero_heat_flux!(model)
end

end # module LithosphericExtensionDepthDependentBC 