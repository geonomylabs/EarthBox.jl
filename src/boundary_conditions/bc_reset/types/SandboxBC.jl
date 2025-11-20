module SandboxBC

import EarthBox.ModelDataContainer: ModelData
import ...SetBoundaryConditions

"""
    set_bcs!(model::ModelData)

Boundary condition coefficients for sandbox shortening experiment.

Boundary conditions for sandbox shortening experiment. No slip on bottom and 
left boundaries. Free slip on top and right side. prescribed internal boundary 
conditions and boundary friction.
"""
function set_bcs!(model::ModelData)
    SetBoundaryConditions.Pressure.upper_left_corner!(model)
    SetBoundaryConditions.TopBoundary.free_slip!(model)
    SetBoundaryConditions.BottomBoundary.no_slip!(model)
    SetBoundaryConditions.LeftBoundary.no_slip!(model)
    SetBoundaryConditions.RightBoundary.free_slip!(model)
    SetBoundaryConditions.Internal.mobile_wall_sandbox!(model)
    SetBoundaryConditions.TopBoundary.temperature!(model)
    SetBoundaryConditions.BottomBoundary.temperature!(model)
    SetBoundaryConditions.LeftBoundary.zero_heat_flux!(model)
    SetBoundaryConditions.RightBoundary.zero_heat_flux!(model)
end

end # module SandboxBC 