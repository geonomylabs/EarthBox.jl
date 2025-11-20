module VerticalChannelFlowShearHeatingBC

import EarthBox.ModelDataContainer: ModelData
import ...SetBoundaryConditions: OpenVerticalChannel
import ...SetBoundaryConditions: LeftBoundary
import ...SetBoundaryConditions: RightBoundary
import ...SetBoundaryConditions: TopBoundary
import ...SetBoundaryConditions: BottomBoundary

"""
    set_bcs!(model::ModelData)

Set boundary condition coefficients for channel flow with shear heating benchmark.

Vertical channel flow with prescribed pressure and velocity gradients
along upper and lower boundaries. No-slip bc's on sides. Used for
channel flow benchmark with variable thermal conductivity and shear heating.

Note this is identical to bc14. Consolidate options.
"""
function set_bcs!(model::ModelData)
    OpenVerticalChannel.set_open_vertical_channel!(model)
    LeftBoundary.no_slip!(model)
    RightBoundary.no_slip!(model)
    TopBoundary.constant_vertical_temperature_gradient!(model, dtdy=0.0)
    BottomBoundary.constant_vertical_temperature_gradient!(model, dtdy=0.0)
    LeftBoundary.variable_temperature!(model)
    RightBoundary.variable_temperature!(model)
end

end # module VerticalChannelFlowShearHeatingBC 