module LavaFlowSolverManager

include("core/MakeFlow.jl")

import EarthBox.PrintFuncs: print_flow_info
import EarthBox.ModelDataContainer: ModelData
import EarthBox.MathTools: linear_interp_single, generate_normal_random_number
import EarthBox.ModelStructureManager.TopAndBottom: calculate_top_and_bottom_of_swarm_opt
import EarthBox.DataStructures: SedimentTransportParameters
import EarthBox.EBCopy: copy_array_1d!
import EarthBox.Compaction.ApplyCompaction: apply_compaction_model
import .MakeFlow: make_flow
import ..PrintLavaFlowInfo: print_flow_info

mutable struct LavaFlowSolver
    topo_gridx::Vector{Float64}
    topo_gridy::Vector{Float64}
    topo_gridy_initial::Vector{Float64}
    sediment_and_flow_thickness_initial::Union{Vector{Float64}, Nothing}
    sediment_and_flow_thickness_initial_compacted::Union{Vector{Float64}, Nothing}
    sediment_and_flow_thickness_total::Union{Vector{Float64}, Nothing}
    compaction_displacement_max::Union{Vector{Float64}, Nothing}
    total_lava_thickness::Union{Vector{Float64}, Nothing}
    eruption_location_x_min::Float64
    width_eruption_domain::Float64
    total_extrusion_volume::Float64
    number_of_flows_per_model_time_step::Int
    residual_lava_thickness_subaerial::Float64
    residual_laval_thickness_submarine::Float64
    y_sealevel::Float64
    lava_flow_decompaction_parameters::SedimentTransportParameters
    use_random_eruption_location::Bool
    use_normal_eruption_location::Bool
    use_compaction_correction::Bool
    decimation_factor::Int
end

function LavaFlowSolver(
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64},
    sediment_and_flow_thickness_initial::Union{Vector{Float64}, Nothing},
    eruption_location_x_min::Float64,
    width_eruption_domain::Float64,
    total_extrusion_volume::Float64,
    number_of_flows_per_model_time_step::Int,
    residual_lava_thickness_subaerial::Float64,
    residual_laval_thickness_submarine::Float64,
    y_sealevel::Float64,
    lava_flow_decompaction_parameters::SedimentTransportParameters;
    use_random_eruption_location::Bool=false,
    use_normal_eruption_location::Bool=false,
    use_compaction_correction::Bool=false,
    decimation_factor::Int=1
)
    LavaFlowSolver(
        topo_gridx,
        topo_gridy,
        copy(topo_gridy),
        sediment_and_flow_thickness_initial,
        nothing,
        nothing,
        nothing,
        nothing,
        eruption_location_x_min,
        width_eruption_domain,
        total_extrusion_volume,
        number_of_flows_per_model_time_step,
        residual_lava_thickness_subaerial,
        residual_laval_thickness_submarine,
        y_sealevel,
        lava_flow_decompaction_parameters,
        use_random_eruption_location,
        use_normal_eruption_location,
        use_compaction_correction,
        decimation_factor
    )
end

""" Extrude magma at the surface and apply compaction correction.

    Updates the y-coordinates of the topography grid (meters).
"""
function extrude_magma(
    solver::LavaFlowSolver,
    model::Union{ModelData, Nothing}=nothing
)::Nothing
    print_extrusion_info(solver)

    (
        total_lava_thickness, 
        sediment_and_flow_thickness_total,
        sediment_and_flow_thickness_initial_compacted
    ) = extrude_magma_at_surface(
        model, 
        solver.topo_gridx, 
        solver.topo_gridy,
        solver.sediment_and_flow_thickness_initial,
        solver.eruption_location_x_min, 
        solver.width_eruption_domain,
        solver.total_extrusion_volume,
        solver.number_of_flows_per_model_time_step,
        solver.residual_lava_thickness_subaerial,
        solver.residual_laval_thickness_submarine,
        solver.y_sealevel, solver.use_random_eruption_location,
        solver.use_normal_eruption_location, solver.decimation_factor,
        solver.lava_flow_decompaction_parameters,
        solver.use_compaction_correction
    )

    solver.total_lava_thickness = total_lava_thickness
    solver.sediment_and_flow_thickness_total = sediment_and_flow_thickness_total
    solver.sediment_and_flow_thickness_initial_compacted = sediment_and_flow_thickness_initial_compacted
    return nothing
end

function print_extrusion_info(solver::LavaFlowSolver)
    print_flow_info("Extruding magma at the surface.")
    print_flow_info("Number of flows per model time step: $(solver.number_of_flows_per_model_time_step)")
end

""" Extrude magma at the surface.

# Updated Arrays
- topo_gridy: Vector{Float64}
    The y-coordinates of the topography grid (meters).

# Returns
- total_lava_thickness: Vector{Float64} (xnum)
    The total thickness of new lava on the topography grid (meters). This
    thickness represent the total new thickness of lava extruded at the
    surface for this time step. The total lava thickness is used to update
    the extrusion thickness array in the topography grid array
    (gridt[6, i]), which inturn is used to transform markers to volcanic
    flow markers.

- sediment_and_flow_thickness_total: Vector{Float64} (xnum)
    The total thickness of sediment and decompacted lava on the topography
    grid (meters). This includes compacted pre-existing sediment and flows
    and decompacted new flows.

- sediment_and_flow_thickness_initial_compacted: Vector{Float64} (xnum)
    The initial sediment and volcanics thickness after compaction (meters).
    This does not include the new flow thickness.

"""
function extrude_magma_at_surface(
    model::Union{ModelData, Nothing},
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64},
    sediment_and_flow_thickness_initial::Vector{Float64},
    eruption_location_x_min::Float64,
    width_eruption_domain::Float64,
    total_extrusion_volume::Float64,
    number_of_flows_per_model_time_step::Int,
    residual_lava_thickness_subaerial::Float64,
    residual_laval_thickness_submarine::Float64,
    y_sealevel::Float64,
    use_random_eruption_location::Bool,
    use_normal_eruption_location::Bool,
    decimation_factor::Int,
    lava_flow_decompaction_parameters::SedimentTransportParameters,
    use_compaction_correction::Bool=false
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
    topo_gridy_initial = copy(topo_gridy)

    (
        eruption_location_out_of_bounds, 
        total_lava_thickness
    ) = lava_flow_loop(
        topo_gridx, topo_gridy, eruption_location_x_min, width_eruption_domain,
        total_extrusion_volume, number_of_flows_per_model_time_step,
        residual_lava_thickness_subaerial, residual_laval_thickness_submarine,
        y_sealevel, use_random_eruption_location, use_normal_eruption_location,
        decimation_factor
    )

    sediment_and_flow_thickness_total = copy(sediment_and_flow_thickness_initial)
    sediment_and_flow_thickness_initial_compacted = copy(sediment_and_flow_thickness_initial)

    if use_compaction_correction && !eruption_location_out_of_bounds && model !== nothing
        (
            total_lava_thickness, 
            sediment_and_flow_thickness_total,
            sediment_and_flow_thickness_initial_compacted
         ) = apply_compaction_model(
            model, topo_gridx, topo_gridy, topo_gridy_initial,
            sediment_and_flow_thickness_initial,
            lava_flow_decompaction_parameters
        )
    end

    return (
        total_lava_thickness, 
        sediment_and_flow_thickness_total,
        sediment_and_flow_thickness_initial_compacted
        )
end

""" Extrude all lava flows for this time step.

# Arguments
- topo_gridx: Vector{Float64} (xnum)
    The x-coordinates of the topography grid (meters).

- topo_gridy: Vector{Float64} (xnum)
    The y-coordinates of the topography grid (meters).

- eruption_location_x_min: Float64
    The minimum eruption location in x-direction (meters).

- width_eruption_domain: Float64
    The width of the eruption domain (meters).

- total_extrusion_volume: Float64
    The total extrusion volume (m^3).

- number_of_flows_per_model_time_step: Int
    The number of flows per model time step.

- residual_lava_thickness_subaerial: Float64
    The residual lava thickness for subaerial flows (meters).

- residual_lava_thickness_submarine: Float64
    The residual lava thickness for submarine flows (meters).

- y_sealevel: Float64
    The y-coordinate of the sea level (meters).

- use_random_eruption_location: Bool
    A boolean flag to use a random eruption location.

- use_normal_eruption_location: Bool
    A boolean flag to use a normal distribution for the eruption location.

- decimation_factor: Int
    The decimation factor for the topography grid.


# Updated Array
- topo_gridy: Array((xnum), dtype=np.float64)
    The y-coordinates of the topography grid (meters).

# Returns
- eruption_location_out_of_bounds: Bool

- total_lava_thickness: Array((xnum), dtype=np.float64)
    The total thickness of new lava on the topography grid (meters). This
    thickness represent the total new thickness of lava extruded at the
    surface for this time step. The total lava thickness is used to update
    the extrusion thickness array in the topography grid array
    (gridt[6, i]), which inturn is used to transform markers to volcanic
    flow markers.

"""
function lava_flow_loop(
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64},
    eruption_location_x_min::Float64,
    width_eruption_domain::Float64,
    total_extrusion_volume::Float64,
    number_of_flows_per_model_time_step::Int,
    residual_lava_thickness_subaerial::Float64,
    residual_lava_thickness_submarine::Float64,
    y_sealevel::Float64,
    use_random_eruption_location::Bool=true,
    use_normal_eruption_location::Bool=false,
    decimation_factor::Int=1
)::Tuple{Bool, Vector{Float64}}
    xnum = length(topo_gridx)
    xmin = topo_gridx[1]
    xmax = topo_gridx[end]

    flow_volume = total_extrusion_volume / number_of_flows_per_model_time_step

    total_lava_thickness_compacted = zeros(Float64, xnum)
    eruption_location_out_of_bounds = false

    for mm in 1:number_of_flows_per_model_time_step
        (
            eruption_location_x, 
            eruption_location_out_of_bounds
        ) = calculate_eruption_x_location(
                eruption_location_x_min, width_eruption_domain,
                use_random_eruption_location, use_normal_eruption_location,
                eruption_location_out_of_bounds, topo_gridx
            )

        (
            eruption_stype, 
            residual_lava_thickness, 
            eruption_location_y
        ) = determine_eruption_style(
                topo_gridx, topo_gridy, eruption_location_x,
                y_sealevel, residual_lava_thickness_subaerial,
                residual_lava_thickness_submarine
            )

        flow_thickness = zeros(Float64, xnum)
        if xmin < eruption_location_x < xmax
            make_flow(
                topo_gridx, topo_gridy, flow_thickness, flow_volume,
                residual_lava_thickness, eruption_location_x;
                decimation_factor=decimation_factor, tolerance=1e-4, nmax=1000,
                use_single_pulse=false
            )
        end

        update_topography_with_flow_thickness(topo_gridy, flow_thickness)
        update_total_lava_thickness_compacted(
            total_lava_thickness_compacted, flow_thickness)

        print_info = false
        if print_info
            print_flow_info(
                mm, eruption_location_x, eruption_location_y, y_sealevel,
                eruption_stype, residual_lava_thickness, flow_volume,
                maximum(flow_thickness)
            )
        end
    end

    total_lava_thickness = copy(total_lava_thickness_compacted)
    return eruption_location_out_of_bounds, total_lava_thickness
end

function calculate_eruption_x_location(
    eruption_location_x_min::Float64,
    width_eruption_domain::Float64,
    use_random_eruption_location::Bool,
    use_normal_eruption_location::Bool,
    eruption_location_out_of_bounds::Bool,
    topo_gridx::Vector{Float64}
)
    if eruption_location_x_min < topo_gridx[1]
        eruption_location_out_of_bounds = true
        eruption_location_x = eruption_location_x_min
    else
        eruption_location_x = calculate_eruption_location(
            eruption_location_x_min, width_eruption_domain,
            use_random_eruption_location=use_random_eruption_location, 
            use_normal_eruption_location=use_normal_eruption_location
        )
    end
    return eruption_location_x, eruption_location_out_of_bounds
end

""" Determine the eruption style.

# Returns
- eruption_stype: String
    The eruption style. This can be either 'subaerial' or 'submarine'.
- residual_lava_thickness: Float64
    The residual lava thickness (meters) based on style.
- eruption_location_y: Float64
    The y-coordinate of the eruption location.
"""
function determine_eruption_style(
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64},
    eruption_location_x::Float64,
    y_sealevel::Float64,
    residual_lava_thickness_subaerial::Float64,
    residual_lava_thickness_submarine::Float64
)::Tuple{String, Float64, Float64}
    eruption_location_y = linear_interp_single(
        topo_gridx, topo_gridy, eruption_location_x)
    if eruption_location_y > y_sealevel
        residual_lava_thickness = residual_lava_thickness_submarine
        eruption_stype = "submarine"
    else
        residual_lava_thickness = residual_lava_thickness_subaerial
        eruption_stype = "subaerial"
    end
    return eruption_stype, residual_lava_thickness, eruption_location_y
end

function calculate_eruption_location(
    eruption_location_x_min::Float64,
    width_eruption_domain::Float64;
    use_random_eruption_location::Bool=true,
    use_normal_eruption_location::Bool=false
)::Float64
    if use_normal_eruption_location
        eruption_location_x = calculate_probabilistic_eruption_location(
            eruption_location_x_min, width_eruption_domain)
    elseif use_random_eruption_location
        eruption_location_x = eruption_location_x_min + rand() * width_eruption_domain
    else
        eruption_location_x = eruption_location_x_min + width_eruption_domain / 2.0
    end
    return eruption_location_x
end

""" Calculate the eruption location using a normal distribution.
"""
function calculate_probabilistic_eruption_location(
    eruption_location_x_min::Float64,
    width_eruption_domain::Float64
)
    mean = eruption_location_x_min + width_eruption_domain / 2.0
    std_dev = width_eruption_domain / 4.0
    eruption_location_x = generate_normal_random_number(
        mean=mean, standard_deviation=std_dev)
    while eruption_location_x < eruption_location_x_min ||
          eruption_location_x > eruption_location_x_min + width_eruption_domain
        eruption_location_x = generate_normal_random_number(
            mean=mean, standard_deviation=std_dev)
    end
    return eruption_location_x
end

function update_topography_with_flow_thickness(
    topo_gridy::Vector{Float64},
    flow_thickness::Vector{Float64}
)::Nothing
    xnum = length(topo_gridy)
    for i in 1:xnum
        topo_gridy[i] -= flow_thickness[i]
    end
    return nothing
end

""" Update total compacted lava thickness with the new lava flow thickness.
"""
function update_total_lava_thickness_compacted(
    total_lava_thickness_compacted::Vector{Float64},
    flow_thickness::Vector{Float64}
)::Nothing
    xnum = length(total_lava_thickness_compacted)
    for i in 1:xnum
        total_lava_thickness_compacted[i] += flow_thickness[i]
    end
    return nothing
end

""" Update total lava thickness with the new lava flow thickness.

# Updated Array
- total_lava_thickness: Vector{Float64}
    The total decompacted thickness of lava on the topography grid (meters).
"""
function update_total_lava_thickness(
    total_lava_thickness::Vector{Float64},
    flow_thickness_decompacted::Vector{Float64}
)::Nothing
    xnum = length(total_lava_thickness)
    xnum2 = length(flow_thickness_decompacted)
    if xnum != xnum2
        error("The size of the arrays do not match")
    end
    for i in 1:xnum
        total_lava_thickness[i] = (
            total_lava_thickness[i] + flow_thickness_decompacted[i]
        )
    end
    return nothing
end

end # module 