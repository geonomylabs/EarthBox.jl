module ElasticSlabBC

import EarthBox.ModelDataContainer: ModelData
import ...SetBoundaryConditions

"""
    set_bcs!(model::ModelData)

Boundary condition coefficients for elastic slab benchmark.

No slip boundary conditions are applied on the left boundary. Free slip boundary 
conditions are applied on top, bottom and right boundaries.
"""
function set_bcs!(model::ModelData)
    SetBoundaryConditions.Pressure.upper_left_corner!(model)
    SetBoundaryConditions.TopBoundary.free_slip!(model)
    SetBoundaryConditions.BottomBoundary.free_slip!(model)
    SetBoundaryConditions.LeftBoundary.no_slip!(model)
    SetBoundaryConditions.RightBoundary.free_slip!(model)
    SetBoundaryConditions.TopBoundary.temperature!(model)
    SetBoundaryConditions.BottomBoundary.temperature!(model)
    SetBoundaryConditions.LeftBoundary.zero_heat_flux!(model)
    SetBoundaryConditions.RightBoundary.zero_heat_flux!(model)
end

end # module ElasticSlab 