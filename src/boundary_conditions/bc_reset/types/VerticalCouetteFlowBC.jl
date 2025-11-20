module VerticalCouetteFlowBC

import EarthBox.ModelDataContainer: ModelData
import ...SetBoundaryConditions: OpenVerticalChannel
import ...SetBoundaryConditions: LeftBoundary
import ...SetBoundaryConditions: RightBoundary
import ...SetBoundaryConditions: TopBoundary
import ...SetBoundaryConditions: BottomBoundary

"""
    set_bcs!(model::ModelData)

Set boundary condition coefficients for vertical Couette flow (simple shear).

Boundary conditions for vertical Couette flow (simple shear). Open channel
boundary conditions are applied along the top and bottom boundaries to
model vertical flow in an infinite channel. This is accomplished by
prescribing zero vertical velocity gradients along top and bottom
boundaries and by prescribing a vertical pressure gradient, which for this
case is zero. No-slip bc's on are applied on side boundaries with
prescribed y-velocity along the right boundary to induce shear flow.
Zero vertical temperature gradients are prescribed along top and bottom
boundaries. Prescribed temperature is applied along left boundary. An
insulating boundary condition (zero lateral heat flux) is implemented along
the right boundary.
"""
function set_bcs!(model::ModelData)
    OpenVerticalChannel.set_open_vertical_channel!(model)
    LeftBoundary.no_slip!(model)
    RightBoundary.no_slip_with_shear!(model)
    TopBoundary.constant_vertical_temperature_gradient!(model, dtdy=0.0)
    BottomBoundary.constant_vertical_temperature_gradient!(model, dtdy=0.0)
    LeftBoundary.temperature!(model)
    RightBoundary.zero_heat_flux!(model)
end

end # module VerticalCouetteFlowBC 