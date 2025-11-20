module VerticalChannelFlowIsoviscousBC

import EarthBox.ModelDataContainer: ModelData
import ...SetBoundaryConditions: OpenVerticalChannel
import ...SetBoundaryConditions: LeftBoundary
import ...SetBoundaryConditions: RightBoundary
import ...SetBoundaryConditions: TopBoundary
import ...SetBoundaryConditions: BottomBoundary

"""
    set_bcs!(model::ModelData)

Set boundary condition coefficients for isoviscous channel flow benchmark.

Vertical channel flow with prescribed pressure and velocity gradients
along upper and lower boundaries. No-slip bc's on sides.
Prescribed temperatures along top and bottom boundaries. Prescribed
gradients along side boundaries. Special regridding for channel flow.

Used for isoviscous channel flow with thermal advection benchmark with
temperature prescribed along side boundaries.
"""
function set_bcs!(model::ModelData)
    OpenVerticalChannel.set_open_vertical_channel!(model)
    LeftBoundary.no_slip!(model)
    RightBoundary.no_slip!(model)
    TopBoundary.constant_vertical_temperature_gradient!(model)
    BottomBoundary.constant_vertical_temperature_gradient!(model)
    LeftBoundary.variable_temperature!(model)
    RightBoundary.variable_temperature!(model)
end

end # module VerticalChannelFlowIsoviscousBC 