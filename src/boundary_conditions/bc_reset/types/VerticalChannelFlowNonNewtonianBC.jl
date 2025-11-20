module VerticalChannelFlowNonNewtonianBC

import EarthBox.ModelDataContainer: ModelData
import ...SetBoundaryConditions: OpenVerticalChannel
import ...SetBoundaryConditions: LeftBoundary
import ...SetBoundaryConditions: RightBoundary
import ...SetBoundaryConditions: TopBoundary
import ...SetBoundaryConditions: BottomBoundary

"""
    set_bcs!(model::ModelData)

Set boundary condition coefficients for Non-Newtonian channel flow benchmark.

Vertical channel flow with open channel boundary condition with 
prescribed pressure and velocity gradients along upper and lower
boundaries. No-slip bc's on sides. Prescribed temperatures along top and 
bottom boundaries. Prescribed gradients along side boundaries. Special 
regridding for channel flow.

This set of boundary conditions is used for Non-Newtonian channel flow 
benchmark.
"""
function set_bcs!(model::ModelData)
    OpenVerticalChannel.set_open_vertical_channel!(model)
    LeftBoundary.no_slip!(model)
    RightBoundary.no_slip!(model)
    TopBoundary.temperature!(model)
    BottomBoundary.temperature!(model)
    LeftBoundary.zero_heat_flux!(model)
    RightBoundary.zero_heat_flux!(model)
end

end # module VerticalChannelFlowNonNewtonianBC 