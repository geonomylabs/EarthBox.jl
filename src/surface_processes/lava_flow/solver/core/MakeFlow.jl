""" Make lava flow by injecting flow pulses
"""
module MakeFlow

include("LavaFlowPulse.jl")

import .LavaFlowPulse: lava_flow_pulse
import .LavaFlowPulse: radiate_indices!
import EarthBox.MathTools: linear_interp_vals!

""" Extrude lava flow on to topography.

Convenience wrapper that allocates the 6 decimation/pulse buffers
needed by `make_flow!` and forwards. Use `make_flow!` directly when
called from a tight loop (e.g. `lava_flow_loop`) so the buffers can
be preallocated by the caller and reused across flows.

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
    indices = 1:decimation_factor:length(topo_gridx_orig)
    xnum_decimated = length(indices)
    topo_gridx_decimated     = Vector{Float64}(undef, xnum_decimated)
    topo_gridy_decimated     = Vector{Float64}(undef, xnum_decimated)
    flow_thickness_decimated = Vector{Float64}(undef, xnum_decimated)
    lava_thickness_old       = Vector{Float64}(undef, xnum_decimated)
    sorted_indices           = Vector{Int}(undef,     xnum_decimated)
    distances_scratch        = Vector{Int}(undef,     xnum_decimated)
    make_flow!(
        topo_gridx_orig, topo_gridy_orig, flow_thickness_orig,
        topo_gridx_decimated, topo_gridy_decimated, flow_thickness_decimated,
        lava_thickness_old, sorted_indices, distances_scratch,
        total_flow_volume, residual_lava_thickness, x_eruption_location;
        decimation_factor=decimation_factor,
        tolerance=tolerance,
        nmax=nmax,
        use_single_pulse=use_single_pulse,
    )
    return nothing
end

""" Allocation-free variant of `make_flow`. The 6 buffer arguments are
filled and reused by this function; their sizes must equal
`length(1:decimation_factor:length(topo_gridx_orig))` (the decimated
grid size). Caller is responsible for sizing them correctly.

The buffers are:
- `topo_gridx_decimated`, `topo_gridy_decimated`, `flow_thickness_decimated`:
    Float64 vectors that receive the strided copy from the corresponding
    `_orig` arrays. Mutated each call.
- `lava_thickness_old`: Float64 scratch buffer for the convergence loop
    inside `lava_flow_pulse`. Contents on entry irrelevant.
- `sorted_indices`, `distances_scratch`: Int buffers used by
    `radiate_indices!` to compute the radiate-from-eruption permutation.
"""
function make_flow!(
    topo_gridx_orig::Vector{Float64},
    topo_gridy_orig::Vector{Float64},
    flow_thickness_orig::Vector{Float64},
    topo_gridx_decimated::Vector{Float64},
    topo_gridy_decimated::Vector{Float64},
    flow_thickness_decimated::Vector{Float64},
    lava_thickness_old::Vector{Float64},
    sorted_indices::Vector{Int},
    distances_scratch::Vector{Int},
    total_flow_volume::Float64,
    residual_lava_thickness::Float64,
    x_eruption_location::Float64;
    decimation_factor::Int=1,
    tolerance::Float64=1e-4,
    nmax::Int=1000,
    use_single_pulse::Bool=false
)::Nothing
    # Decimate the grid via explicit copy into the caller-provided
    # buffers — eliminates the 3 fancy-index-slice allocations that the
    # original `topo_gridx_orig[indices]` formulation incurred.
    indices = 1:decimation_factor:length(topo_gridx_orig)
    @inbounds for (k, idx) in enumerate(indices)
        topo_gridx_decimated[k]     = topo_gridx_orig[idx]
        topo_gridy_decimated[k]     = topo_gridy_orig[idx]
        flow_thickness_decimated[k] = flow_thickness_orig[idx]
    end

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
        # `eruption_node` is invariant within a flow (x_eruption_location
        # is fixed), so `radiate_indices!` is called once and
        # `sorted_indices` is reused for every pulse.
        xnum_decimated = length(topo_gridx_decimated)
        eruption_node = floor(Int, x_eruption_location / dx) + 1
        if eruption_node < 1 || eruption_node > xnum_decimated
            eruption_node = 1
        end
        radiate_indices!(
            sorted_indices, distances_scratch, xnum_decimated, eruption_node)
        for _ in 1:npulses
            lava_flow_pulse(
                topo_gridx_decimated,
                topo_gridy_decimated,
                flow_thickness_decimated,
                lava_thickness_old,
                sorted_indices,
                eruption_node,
                pulse_thickness,
                residual_lava_thickness;
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