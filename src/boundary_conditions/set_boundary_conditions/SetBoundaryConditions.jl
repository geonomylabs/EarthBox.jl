module SetBoundaryConditions

include("utils/SetBCTools.jl")
include("heat_boundary_conditions/ConstantGradientBC.jl")
include("heat_boundary_conditions/TemperatureBC.jl")
include("heat_boundary_conditions/ZeroHeatFluxBC.jl")
include("stokes_boundary_conditions/PressureBC.jl")
include("stokes_boundary_conditions/StokesGradientBC.jl")
include("stokes_boundary_conditions/VelocityBC.jl")
include("stokes_boundary_conditions/VelocityInflowOutflow.jl")
include("stokes_boundary_conditions/VerticalChannelBC.jl")
include("stokes_boundary_conditions/NoSlipBC.jl")
include("stokes_boundary_conditions/FreeSlipBC.jl")
include("stokes_boundary_conditions/InflowOutflowBC.jl")
include("stokes_boundary_conditions/InternalBC.jl")
include("domains/TopBoundary.jl")
include("domains/BottomBoundary.jl")
include("domains/LeftBoundary.jl")
include("domains/RightBoundary.jl")
include("domains/Internal.jl")
include("domains/Pressure.jl")
include("domains/OpenVerticalChannel.jl")

end # module SetBoundaryConditions 