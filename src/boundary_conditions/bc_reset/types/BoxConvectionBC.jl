module BoxConvectionBC

import EarthBox.ModelDataContainer: ModelData
import ...SetBoundaryConditions

"""
    set_bcs!(model::ModelData)

Boundary condition coefficients for convection benchmark.

Boundary conditions for convection in a box. Free slip on all sides with 
prescribed temperature on top and bottom, and insulating bc's on sides.
"""
function set_bcs!(model::ModelData)
    SetBoundaryConditions.Pressure.upper_left_corner!(model)
    SetBoundaryConditions.TopBoundary.free_slip!(model)
    SetBoundaryConditions.BottomBoundary.free_slip!(model)
    SetBoundaryConditions.LeftBoundary.free_slip!(model)
    SetBoundaryConditions.RightBoundary.free_slip!(model)
    SetBoundaryConditions.TopBoundary.temperature!(model)
    SetBoundaryConditions.BottomBoundary.temperature!(model)
    SetBoundaryConditions.LeftBoundary.zero_heat_flux!(model)
    SetBoundaryConditions.RightBoundary.zero_heat_flux!(model)
end

end # module