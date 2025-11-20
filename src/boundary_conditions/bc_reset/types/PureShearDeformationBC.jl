module PureShearDeformationBC

import EarthBox.ModelDataContainer: ModelData
import ...SetBoundaryConditions

"""
    set_bcs!(model::ModelData)

Boundary condition coefficients for pure shear benchmark.

Boundary conditions for pure shear experiment. Free slip with orthogonal 
prescribed velocity. Temperature prescribed on top, bottom and side boundaries. 
No marker recycling.
"""
function set_bcs!(model::ModelData)
    # Unpack model data
    vx = model.bcs.parameters.velocity.full_velocity_extension.value/2.0
    vy = model.bcs.parameters.velocity.full_velocity_extension.value/2.0
    
    SetBoundaryConditions.Pressure.upper_left_corner!(model)
    SetBoundaryConditions.TopBoundary.inflow_outflow!(model, vy)
    SetBoundaryConditions.BottomBoundary.inflow_outflow!(model, -vy)
    SetBoundaryConditions.LeftBoundary.inflow_outflow!(model, -vx)
    SetBoundaryConditions.RightBoundary.inflow_outflow!(model, vx)
    SetBoundaryConditions.TopBoundary.temperature!(model)
    SetBoundaryConditions.BottomBoundary.temperature!(model)
    SetBoundaryConditions.LeftBoundary.temperature!(model)
    SetBoundaryConditions.RightBoundary.temperature!(model)
end

end # module PureShearDeformation 