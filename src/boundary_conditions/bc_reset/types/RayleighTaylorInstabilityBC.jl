module RayleighTaylorInstabilityBC

import EarthBox.ModelDataContainer: ModelData
import ...SetBoundaryConditions

"""
    set_bcs!(model::ModelData)

Boundary condition coefficients for Rayleigh-Taylor instability.

Boundary conditions for Ramberg benchmark. No slip on top and bottom and free 
slip on sides.
"""
function set_bcs!(model::ModelData)
    SetBoundaryConditions.Pressure.upper_left_corner!(model)
    SetBoundaryConditions.TopBoundary.no_slip!(model)
    SetBoundaryConditions.BottomBoundary.no_slip!(model)
    SetBoundaryConditions.LeftBoundary.free_slip!(model)
    SetBoundaryConditions.RightBoundary.free_slip!(model)
    SetBoundaryConditions.TopBoundary.temperature!(model)
    SetBoundaryConditions.BottomBoundary.temperature!(model)
    SetBoundaryConditions.LeftBoundary.temperature!(model)
    SetBoundaryConditions.RightBoundary.temperature!(model)
end

end # module RayleighTaylorInstabilityBC 