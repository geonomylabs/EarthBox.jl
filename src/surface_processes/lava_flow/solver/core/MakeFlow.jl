""" Make lava flow by injecting flow pulses
"""
module MakeFlow

include("LavaFlowPulse.jl")

import .LavaFlowPulse: lava_flow_pulse
import EarthBox.MathTools: linear_interp_vals!

""" Extrude lava flow on to topography.

# Arguments
- `topo_gridx_orig`: The x-coordinates (meters) of the topography grid
- `topo_gridy_orig`: The y-coordinates (meters) of the topography grid
- `flow_thickness_orig`: The thickness of lava on the topography grid (meters)
- `total_flow_volume`: Total volume of lava flow
- `residual_lava_thickness`: The residual thickness of lava (meters)
- `x_eruption_location`: The x-coordinate of the eruption point (meters)
- `decimation_factor`: Factor to decimate the grid
- `tolerance`: The tolerance for the lava flow model
- `nmax`: The maximum number of iterations for the lava flow model
- `use_single_pulse`: Whether to use a single pulse

"""
function make_flow(
    topo_gridx_orig::Vector{Float64},
    topo_gridy_orig::Vector{Float64},
    flow_thickness_orig::Vector{Float64},
    total_flow_volume::Float64,
    residual_lava_thickness::Float64,
    x_eruption_location::Float64;
    decimation_factor::Int=1,
    tolerance::Float64=1e-4,
    nmax::Int=1000,
    use_single_pulse::Bool=false
)::Nothing
    # Decimate the grid
    indices = 1:decimation_factor:length(topo_gridx_orig)
    topo_gridx_decimated = topo_gridx_orig[indices]
    topo_gridy_decimated = topo_gridy_orig[indices]
    flow_thickness_decimated = flow_thickness_orig[indices]

    dx = topo_gridx_decimated[2] - topo_gridx_decimated[1]
    npulses = calculate_number_of_pulses(
        dx, total_flow_volume, residual_lava_thickness,
        use_single_pulse=use_single_pulse
    )
    pulse_thickness = total_flow_volume / dx / npulses

    print_info = false
    if print_info
        print_flow_info(npulses, pulse_thickness)
    end

    xmin = topo_gridx_decimated[1]
    xmax = topo_gridx_decimated[end]
    if xmin < x_eruption_location < xmax
        for _ in 1:npulses
            _niterations = lava_flow_pulse(
                topo_gridx_decimated,
                topo_gridy_decimated,
                flow_thickness_decimated,
                pulse_thickness,
                residual_lava_thickness,
                x_eruption_location;
                tolerance=tolerance,
                nmax=nmax
            )
        end
    end
    linear_interp_vals!(
        topo_gridx_decimated, flow_thickness_decimated,
        topo_gridx_orig, flow_thickness_orig
    )
    return nothing
end

function print_flow_info(
    npulses::Int,
    pulse_thickness::Float64
)::Nothing
    println(">> Number of pulses for flow event: ", npulses)
    println(">> Thickness per pulse (meters): ", pulse_thickness)
    return nothing
end

""" Calculate the number of pulses for the lava flow model.

# Arguments
- `dx`: Grid spacing
- `total_flow_volume`: Total volume of lava flow
- `residual_lava_thickness`: The residual thickness of lava (meters)
- `use_single_pulse`: Whether to use a single pulse

# Returns
- Number of pulses
"""
function calculate_number_of_pulses(
    dx::Float64,
    total_flow_volume::Float64,
    residual_lava_thickness::Float64;
    use_single_pulse::Bool=false
)::Int
    if use_single_pulse
        npulses = 1
    else
        npulses = floor(Int, total_flow_volume / dx / residual_lava_thickness)
    end
    return max(1, npulses)
end

end # module 