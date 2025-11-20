module SolidBodyRotationBC

import EarthBox.ModelDataContainer: ModelData
import ...SetBoundaryConditions

"""
    set_bcs!(model::ModelData)

Boundary condition coefficients for solid body rotation.

Solid body rotation benchmark with fixed boundaries (Stokes equations are not 
solved). Prescribed gradients along side boundaries.
"""
function set_bcs!(model::ModelData)
    SetBoundaryConditions.Pressure.top_and_bottom!(model)
    SetBoundaryConditions.TopBoundary.no_slip!(model)
    SetBoundaryConditions.BottomBoundary.no_slip!(model)
    SetBoundaryConditions.LeftBoundary.no_slip!(model)
    SetBoundaryConditions.RightBoundary.no_slip!(model)
    SetBoundaryConditions.TopBoundary.constant_vertical_temperature_gradient!(model, dtdy=0.0)
    SetBoundaryConditions.BottomBoundary.constant_vertical_temperature_gradient!(model, dtdy=0.0)
    SetBoundaryConditions.LeftBoundary.temperature!(model)
    SetBoundaryConditions.RightBoundary.temperature!(model)
end

end # module SolidBodyRotationBC 