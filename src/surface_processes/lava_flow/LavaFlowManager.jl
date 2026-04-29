""" Functions used to model lava flows using the Cellular Automata approach.
"""
module LavaFlowManager

include("utils/GetData.jl")
include("utils/PlotLavaFlow.jl")
include("utils/PrintLavaFlowInfo.jl")
include("solver/LavaFlowSolverManager.jl")

import EarthBox.ModelDataContainer: ModelData
import EarthBox.PrintFuncs: print_info
import EarthBox.MathTools: linear_interp_at_x_location
import EarthBox.EBCopy: copy_new_topography_to_topography_array
import EarthBox.EBCopy: copy_topography_coordinate_arrays
import EarthBox.DataStructures: SedimentTransportParameters
import EarthBox.Compaction: get_boolean_options_for_compaction
import .GetData: get_extrusion_parameters
import .GetData: get_extrusion_location_parameters
import .PlotLavaFlow: plot_lava_thickness
import .PrintLavaFlowInfo: print_lava_model_info
import .LavaFlowSolverManager: LavaFlowSolver
import .LavaFlowSolverManager: extrude_magma
import .LavaFlowSolverManager: reset!
import .LavaFlowSolverManager: reset_timestep!
import .LavaFlowSolverManager: is_compatible

function lava_flow_manager(
    model::ModelData,
    sediment_and_flow_thickness_initial::Union{Array{Float64}, Nothing},
    use_magma_flush::Bool
)::Nothing
    iuse_extrusion = model.melting.parameters.extrusion.iuse_extrusion.value
    timesum = model.timestep.parameters.main_time_loop.timesum.value
    sec_per_myr = model.conversion.parameters.sec_per_Myr.value
    timesum_myr = timesum / sec_per_myr
    time_of_next_eruption_myr = model.melting.parameters.extrusion.time_of_next_eruption_myr.value

    if iuse_extrusion == 1 || use_magma_flush
        if timesum_myr >= time_of_next_eruption_myr
            print_info("The time of the next eruption event is $time_of_next_eruption_myr Myr. Running lava flow model at $timesum_myr Myr.", level=2)
            run_lava_flow_model(model, sediment_and_flow_thickness_initial)
            # Reset extrusion volumes to zero
            zero_out_extrusion_volumes(model)
        else
            print_info("The time of the next eruption event is $time_of_next_eruption_myr Myr. Skipping lava flow model at $timesum_myr Myr.", level=2)
        end
    end
    return nothing
end

function zero_out_extrusion_volumes(model::ModelData)::Nothing
    extrusion_volumes = model.melting.arrays.extraction.extrusion_volumes.array
    nvalues = size(extrusion_volumes, 1)
    for idrainage_basin in 1:nvalues
        extrusion_volumes[idrainage_basin] = 0.0
    end
    return nothing
end

""" Run the lava flow model including compaction of old sediment.

Topography stored in gridt[1, i] is updated to account for the lava flow
and compaction of pre-existing sediment. If use_compaction_correction is
set to True, the volume melt that is extruded from the mantle is
decompacted to account for pore space in the lava flow.

If use_compaction_correction is True, old sediment markers are advected
using the compaction displacement field.

# Arguments
- `model::ModelData`: The model data object.
- `sediment_and_flow_thickness_initial::Vector{Float64}`: 
    The initial sediment and flow thickness.

# Updated Arrays
- `gridt::Array{Float64, 2}`: 
    Topography array. This function updates the y-coordinate of topography
    nodes stored in gridt[i, 1] and the extrusion thickness stored in
    gridt[i, 6].
"""
function run_lava_flow_model(
    model::ModelData,
    sediment_and_flow_thickness_initial::Union{Array{Float64}, Nothing}
)::Nothing
    (
        y_sealevel,
        characteristic_flow_length_subaerial,
        characteristic_flow_length_submarine,
        residual_lava_thickness_subaerial,
        residual_laval_thickness_submarine,
        use_random_eruption_location,
        use_normal_eruption_location,
        decimation_factor
    ) = get_extrusion_parameters(model)

    gridt = model.topography.arrays.gridt.array
    zero_out_extrusion_thickness_for_timestep(gridt)

    (topo_gridx, topo_gridy) = copy_topography_coordinate_arrays(gridt)

    # Translate `nothing` to zeros once here so the persistent solver always
    # holds a real Vector{Float64} (matches the convenience constructor's
    # behavior; required by `reset_timestep!` which uses copy_array_1d!).
    sediment_and_flow_thickness_initial_resolved =
        sediment_and_flow_thickness_initial === nothing ?
            zeros(Float64, length(topo_gridx)) : sediment_and_flow_thickness_initial

    lava_flow_decompaction_parameters = get_lava_flow_decompaction_parameters(
        model.melting.parameters.extrusion.porosity_initial_lava_flow.value,
        model.melting.parameters.extrusion.decay_depth_lava_flow.value
    )

    extrusion_volumes = model.melting.arrays.extraction.extrusion_volumes.array

    use_compaction_correction = get_boolean_options_for_compaction(model)

    # Persistent solver: lazy-initialized on the first call, reused on
    # subsequent timesteps via in-place buffer refresh.
    flow_solver = get_or_init_lava_flow_solver!(
        model,
        topo_gridx,
        topo_gridy,
        sediment_and_flow_thickness_initial_resolved,
        residual_lava_thickness_subaerial,
        residual_laval_thickness_submarine,
        y_sealevel,
        lava_flow_decompaction_parameters,
        use_random_eruption_location,
        use_normal_eruption_location,
        use_compaction_correction,
        decimation_factor
    )

    # The persistent solver owns its topo buffers; the local `topo_gridy`
    # diverges from `flow_solver.topo_gridy` once `extrude_magma` mutates
    # the latter, so we use the solver's buffer for forecasting and write-
    # back across the basin loop.
    ndrainage_basin = model.melting.parameters.extraction.ndrainage_basin.value
    for idrainage_basin in 1:ndrainage_basin
        (
            xmid_molten_domain,
            width_eruption_domain,
            eruption_location_x_min
        ) = get_extrusion_location_parameters(model, idrainage_basin)

        (
            eruption_location_x_forecast,
            eruption_location_y_forecast
        ) = forecast_eruption_location(
            eruption_location_x_min,
            width_eruption_domain,
            flow_solver.topo_gridx,
            flow_solver.topo_gridy
        )

        characteristic_volume_per_flow = calculate_characteristic_volume_per_flow(
            y_sealevel,
            eruption_location_y_forecast,
            characteristic_flow_length_subaerial,
            characteristic_flow_length_submarine,
            residual_lava_thickness_subaerial,
            residual_laval_thickness_submarine
        )

        total_extrusion_volume = extrusion_volumes[idrainage_basin]

        number_of_flows_per_model_time_step = calculate_number_of_flows(
            total_extrusion_volume,
            characteristic_volume_per_flow
        )

        print_info = true
        if print_info
            print_lava_model_info(
                idrainage_basin,
                xmid_molten_domain,
                width_eruption_domain,
                eruption_location_x_min,
                eruption_location_x_forecast,
                eruption_location_y_forecast,
                total_extrusion_volume,
                characteristic_volume_per_flow,
                number_of_flows_per_model_time_step
            )
        end

        if total_extrusion_volume > 0.0
            reset!(
                flow_solver,
                eruption_location_x_min,
                width_eruption_domain,
                total_extrusion_volume,
                number_of_flows_per_model_time_step,
                use_compaction_correction
            )

            extrude_magma(flow_solver, model)

            copy_new_topography_to_topography_array(flow_solver.topo_gridy, gridt)

            update_extrusion_thickness(flow_solver.total_lava_thickness, gridt)

            plot_thickness = false
            if plot_thickness
                plot_lava_thickness(
                    model,
                    flow_solver.topo_gridx,
                    flow_solver.total_lava_thickness,
                    idrainage_basin
                )
            end
        end
    end
    return nothing
end

""" Return a `LavaFlowSolver` whose preallocated buffers persist across
    timesteps.

    On the first call (or when buffers are sized inconsistently with the
    current `topo_gridx` length / `decimation_factor`), constructs a fresh
    solver and stashes it on `model.topography.lava_flow_solver`. On
    subsequent calls, reuses the stashed solver after refreshing its
    per-timestep array buffers via `reset_timestep!` and updating its
    per-timestep scalar configuration in place.

    Per-basin scalars are reset separately by `reset!` inside the basin
    loop in `run_lava_flow_model`.
"""
function get_or_init_lava_flow_solver!(
    model::ModelData,
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64},
    sediment_and_flow_thickness_initial::Vector{Float64},
    residual_lava_thickness_subaerial::Float64,
    residual_laval_thickness_submarine::Float64,
    y_sealevel::Float64,
    lava_flow_decompaction_parameters::SedimentTransportParameters,
    use_random_eruption_location::Bool,
    use_normal_eruption_location::Bool,
    use_compaction_correction::Bool,
    decimation_factor::Int
)::LavaFlowSolver
    cached = model.topography.lava_flow_solver
    toponum = length(topo_gridx)
    if cached isa LavaFlowSolver && is_compatible(cached, toponum, decimation_factor)
        # Refresh per-timestep array state in place.
        reset_timestep!(cached, topo_gridx, topo_gridy, sediment_and_flow_thickness_initial)
        # Refresh per-timestep scalar configuration. (Per-basin scalars are
        # set by `reset!` inside the basin loop; we leave those untouched.)
        cached.residual_lava_thickness_subaerial   = residual_lava_thickness_subaerial
        cached.residual_laval_thickness_submarine  = residual_laval_thickness_submarine
        cached.y_sealevel                          = y_sealevel
        cached.lava_flow_decompaction_parameters   = lava_flow_decompaction_parameters
        cached.use_random_eruption_location        = use_random_eruption_location
        cached.use_normal_eruption_location        = use_normal_eruption_location
        cached.use_compaction_correction           = use_compaction_correction
        return cached
    end

    # First-call path (or buffer-size / decimation_factor changed): build a
    # fresh solver. Per-basin scalars start as placeholders; the basin loop
    # overwrites them via `reset!` before the first `extrude_magma`.
    solver = LavaFlowSolver(
        topo_gridx,
        topo_gridy,
        sediment_and_flow_thickness_initial,
        0.0,
        0.0,
        0.0,
        1,
        residual_lava_thickness_subaerial,
        residual_laval_thickness_submarine,
        y_sealevel,
        lava_flow_decompaction_parameters;
        use_random_eruption_location=use_random_eruption_location,
        use_normal_eruption_location=use_normal_eruption_location,
        use_compaction_correction=use_compaction_correction,
        decimation_factor=decimation_factor
    )
    model.topography.lava_flow_solver = solver
    return solver
end

""" Get sediment transport parameters.

Only porosity_initial and depth_decay_term are used in the
SedimentTransportParameters data structure for the lava flow model.
"""
function get_lava_flow_decompaction_parameters(
    porosity_initial::Float64=0.0,
    decay_depth::Float64=2000.0
)::SedimentTransportParameters
    if decay_depth == 0.0
        decay_depth = 1000.0
    end
    decay_depth_term = 1.0 / decay_depth # 1/m
    return SedimentTransportParameters(
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        porosity_initial,
        decay_depth_term
    )
end

""" Forecast the eruption location.

# Returns
- `eruption_location_x::Float64`: The forecasted x-coordinate of the eruption location
- `eruption_location_y::Float64`: The forecasted y-coordinate of the eruption location
"""
function forecast_eruption_location(
    eruption_location_x_min::Float64,
    width_eruption_domain::Float64,
    topo_gridx::Array{Float64},
    topo_gridy::Array{Float64}
)::Tuple{Float64, Float64}
    eruption_location_x = eruption_location_x_min + width_eruption_domain / 2.0
    eruption_location_y = linear_interp_at_x_location(
        eruption_location_x,
        topo_gridx,
        topo_gridy
    )
    return eruption_location_x, eruption_location_y
end

function calculate_characteristic_volume_per_flow(
    y_sealevel::Float64,
    eruption_location_y_forecast::Float64,
    characteristic_flow_length_subaerial::Float64,
    characteristic_flow_length_submarine::Float64,
    residual_lava_thickness_subaerial::Float64,
    residual_laval_thickness_submarine::Float64
)::Float64
    if eruption_location_y_forecast <= y_sealevel
        characteristic_volume_per_flow = 
            characteristic_flow_length_subaerial * residual_lava_thickness_subaerial
    else
        characteristic_volume_per_flow = 
            characteristic_flow_length_submarine * residual_laval_thickness_submarine
    end
    return characteristic_volume_per_flow
end

function calculate_number_of_flows(
    total_extrusion_volume::Float64,
    characteristic_volume_per_flow::Float64
)::Int
    if total_extrusion_volume < characteristic_volume_per_flow
        number_of_flows_per_model_time_step = 1
    else
        number_of_flows_per_model_time_step = 
            floor(Int, total_extrusion_volume / characteristic_volume_per_flow)
    end
    return number_of_flows_per_model_time_step
end

function update_sediment_thickness_initial(
    sediment_thickness_initial::Array{Float64},
    sediment_thickness_total::Array{Float64}
)::Nothing
    xnum = size(sediment_thickness_initial, 1)
    for i in 1:xnum
        sediment_thickness_initial[i] = sediment_thickness_total[i]
    end
    return nothing
end

""" Set extrusion thickness at topography grid node to zero.

Extrusion thickness is set to zero at the beginning of each time step.

# Updated Arrays
- `gridt::Array{Float64,2}`: shape=(7,toponum)
    Topography grid array. This function updates gridt[7,toponum], which
    stores extrusion thickness (m) at topography node.
"""
function zero_out_extrusion_thickness_for_timestep(gridt::Array{Float64,2})::Nothing
    toponum = size(gridt, 2)
    for i in 1:toponum
        gridt[7, i] = 0.0
    end
    return nothing
end

""" Update extrusion thickness at topography grid node.

# Arguments
- `total_lava_thickness::Array{Float64,1}`: shape=(toponum,)
    Total lava thickness (m) at topography grid node that was emplaced
    during time step.
    
# Updated Arrays
- `gridt::Array{Float64,2}`: shape=(7,toponum)
    Topography grid array. This function updates gridt[7,toponum], which
    stores extrusion thickness (m) at topography node.

"""
function update_extrusion_thickness(
    total_lava_thickness::Array{Float64,1},
    gridt::Array{Float64,2}
)::Nothing
    toponum = size(gridt, 2)
    for i in 1:toponum
        gridt[7, i] = gridt[7, i] + total_lava_thickness[i]
    end
    return nothing
end

end # module 