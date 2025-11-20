module SedimentTransportSolverManager

include("core/Collections.jl")
include("core/BuildSys.jl")
include("core/Diffusivity.jl")
include("core/DownstreamDistance.jl")
include("core/WaterDepth.jl")
include("core/MarkerAdvection.jl")
include("core/Solve.jl")

import EarthBox.PrintFuncs: print_info
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ConversionFuncs: seconds_to_years
import EarthBox.ConversionFuncs: meters_per_seconds_to_mm_per_yr
import EarthBox.Compaction.CompactionTools: calculate_final_thickness_after_burial
import EarthBox.Compaction.MarkerCompaction: calculate_swarm_indices_for_sediment_and_sticky
import EarthBox.Compaction.MarkerCompaction: calculate_x_sorted_swarm_indices
import EarthBox.DataStructures: SedimentTransportParameters
import EarthBox.Compaction.CompactionCorrection: apply_compaction_correction
import EarthBox.Compaction.CompactionCorrection: apply_compaction_correction_for_topography_and_markers
import EarthBox.EBCopy: copy_array_1d!
import .DownstreamDistance: calculate_downstream_distances_for_nodes
import .WaterDepth: calculate_water_depth
import .Diffusivity: make_topo_diffusivity_grid
import .Solve: solve_downhill_diffusion
import .Collections: TransportCollections
import .Collections: update_collections!

mutable struct SedimentTransportSolver
    basic_grid_x_dimensions::Tuple{Float64, Float64}
    topo_gridx::Vector{Float64}
    topo_gridy_initial::Vector{Float64}
    topo_gridy::Vector{Float64}
    sediment_thickness_initial::Vector{Float64}
    sediment_thickness_initial_compacted::Union{Vector{Float64}, Nothing}
    sediment_thickness_total_decompacted::Union{Vector{Float64}, Nothing}
    sediment_transport_parameters::Union{SedimentTransportParameters, Nothing}
    compaction_displacement_max::Union{Vector{Float64}, Nothing}
    y_sealevel::Float64
    pelagic_sedimentation_rate::Float64
    pelagic_sedimentation_rate_reduction_factor::Float64
    pelagic_sedimentation_rate_reduction_time::Float64
    transport_timestep::Float64
    nsteps::Int
    use_collections::Bool
    use_print_debug::Bool
    use_constant_diffusivity::Bool
    use_compaction_correction::Bool
    compaction_correction_type::String
    collections::TransportCollections
end

""" Constructor for transport solver.

# Arguments
- `basic_grid_x_dimensions::Tuple{Float64, Float64}`:
    Basic grid x-dimensions (xmin, xmax) (meters).
- `topo_gridx::Vector{Float64}`:
    Topography x-coordinates (meters).
- `topo_gridy_initial::Vector{Float64}`:
    Initial topography y-coordinates (meters).
- `sediment_thickness_initial::Vector{Float64}`:
    Initial sediment thickness (meters).
- `sediment_transport_parameters::Union{SedimentTransportParameters, Nothing}`:
    This data structure is a namedtuple with the following fields:
    - `subaerial_slope_diffusivity::Float64`:
        Subaerial slope diffusivity (m^2/s). Typical value used in the
        literature is 0.25 m^2/yr (7.9e-9 m^2/s) (e.g. Andres-Martinez
        et al., 2019; Armitage et al., 2015).
    - precipitation_rate : float :
        Precipitation rate (m/s). Used to calculate water flux in
        drainage basins. Typical value used in the literature is
        1 m/yr (3.2e-8 m/s) (e.g. Andres-Martinez et al.; 2019;
        Huffman et al., 2009).
    - subaerial_transport_coefficient : float :
        Subaerial discharge transport coefficient. Used to calculate
        effective subaerial fluvial transport diffusivity that includes
        slope diffusivity, precipitation rate and downstream distances.
        Typical values used in the literature are 1e-4 (low transport) to
        1e-2 (high transport) (e.g. Andres-Martinez et al., 2019; Armitage
        et al., 2015).
    - submarine_slope_diffusivity : float :
        Maximum submarine slope diffusivity (m^2/s). Used in diffusivity
        model that decays exponentially with water depth. Typical value
        used in the literature is 100 m^2/yr (3.2e-9 m^2/s) (e.g.
        Andres-Martinez et al., 2019; Kaufman et al., 1991).
    - submarine_diffusion_decay_depth : float :
        Submarine diffusion decay depth (m). Typical value used in the
        literature is 2000 m (e.g. Andres-Martinez et al., 2019; Kaufman
        et al., 1991).
    - transport_timestep : float :
        Transport timestep (s). Typical value used in the literature is
        1000 years (3.1536e10 s) (e.g. Andres-Martinez et al., 2019).
    - transport_model_duration : float
        Transport model duration (s).
    - porosity_initial : float :
        Initial porosity (fraction) of sediments at sticky-sediment
        interface used to calculate fully compacted pelagic depositional
        thickness and for applying the compaction correction with the
        'constant_property' compaction_correction_type.
    - depth_decay_term : float :
        Depth decay term (1/m) of sediment used to calculate
        fully compacted pelagic depositional thickness and for applying
        the compaction correction with the 'constant_property'
        compaction_correction_type.
- `y_sealevel::Float64`:
    Current model sea level y-coordinate (meters).
- `pelagic_sedimentation_rate::Float64`:
    Decompacted pelagic sedimentation rate (m/s). Typical value used
    in the literature are 0.3 mm/yr (syn-rift) to 0.01 mm/yr
    (post-rift) (e.g. Perez-Gussinye et al., 2020).
- `pelagic_sedimentation_rate_reduction_factor::Float64`:
    The pelagic sedimentation rate is divided by this factor after
    the specified pelagic sedimentation rate reduction time.
- `pelagic_sedimentation_rate_reduction_time::Float64`:
    Time in Myr when the pelagic sedimentation rate is reduced by the
    pelagic sedimentation rate reduction factor.
- `use_collections::Bool`:
    Use transport collections.
- `use_print_debug::Bool`:
    Use print debug.
- `use_constant_diffusivity::Bool`:
    If True constant diffusivity is used equal to the subaerial
    slope diffusivity.
- `use_compaction_correction::Bool`:
    If True compaction correction is applied at each transport time
    step.
- `compaction_correction_type::String`:
    Type of compaction correction:
    - 'constant_property': initial porosity and decay depth are
        constant for all materials in sedimentary basin. For this option
        these parameters are defined from the sediment_transport_parameters
        namedtuple.
    - 'variable_property': initial porosity and decay depth are
        variable for different materials in sedimentary basin. For
        this option these parameters are defined from the model data.

"""
function SedimentTransportSolver(
    basic_grid_x_dimensions::Tuple{Float64, Float64},
    topo_gridx::Vector{Float64},
    topo_gridy_initial::Vector{Float64},
    sediment_thickness_initial::Vector{Float64},
    sediment_transport_parameters::Union{SedimentTransportParameters, Nothing},
    y_sealevel::Float64,
    pelagic_sedimentation_rate::Float64;
    pelagic_sedimentation_rate_reduction_factor::Float64=1.0,
    pelagic_sedimentation_rate_reduction_time::Float64=1000.0,
    use_collections::Bool=true,
    use_print_debug::Bool=false,
    use_constant_diffusivity::Bool=false,
    use_compaction_correction::Bool=false,
    compaction_correction_type::String="constant_property"
)
    nsteps, transport_timestep = update_time_step(sediment_transport_parameters)
    
    return SedimentTransportSolver(
        basic_grid_x_dimensions,
        copy(topo_gridx),
        copy(topo_gridy_initial),
        copy(topo_gridy_initial),
        copy(sediment_thickness_initial),
        nothing,
        nothing,
        sediment_transport_parameters,
        nothing,
        y_sealevel,
        pelagic_sedimentation_rate,
        pelagic_sedimentation_rate_reduction_factor,
        pelagic_sedimentation_rate_reduction_time,
        transport_timestep,
        nsteps,
        use_collections,
        use_print_debug,
        use_constant_diffusivity,
        use_compaction_correction,
        compaction_correction_type,
        TransportCollections(use_collections)
    )
end

function update_time_step(
    parameters::SedimentTransportParameters
)::Tuple{Int, Float64}
    nsteps = max(
        1, floor(Int, parameters.transport_model_duration / parameters.transport_timestep))
    transport_time_step = parameters.transport_model_duration / nsteps
    return nsteps, transport_time_step
end

""" Run sediment transport time steps

# Arguments
- `solver::SedimentTransportSolver`:
    Sediment transport solver.
- `model::Union{ModelData, Nothing}`:
    Model data used for variable property compaction.

# Updated Attributes
- `solver.topo_gridy::Vector{Float64}`:
    Topography at grid nodes (meters).
- `solver.collections::TransportCollections`:
    Transport collections.

"""
function run_sediment_transport_time_steps!(
    solver::SedimentTransportSolver, 
    model::Union{ModelData, Nothing}=nothing
)::Nothing
    (
        downstream_distances_x, drainage_divides_x
    ) = calculate_downstream_distances_for_nodes(solver.topo_gridx, solver.topo_gridy)
    sealevel_x = fill(solver.y_sealevel, length(solver.topo_gridx))
    water_depth_x = calculate_water_depth(solver.topo_gridy, sealevel_x)
    update_collections!(solver.collections, 0, solver.topo_gridy, drainage_divides_x, water_depth_x)
    if solver.compaction_correction_type == "variable_property" && !isnothing(model)
        markers_indices_sedimentary_basin, markers_indices_sticky = 
            calculate_swarm_indices_for_sediment_and_sticky(model)
        x_sorted_marker_indices_sticky = calculate_x_sorted_swarm_indices(
            model, markers_indices_sticky)
        x_sorted_marker_indices_sedimentary_basin = calculate_x_sorted_swarm_indices(
            model, markers_indices_sedimentary_basin)
    end

    pelagic_sedimentation_rate = solver.pelagic_sedimentation_rate
    if !isnothing(model)
        timesum_seconds = model.timestep.parameters.main_time_loop.timesum.value
        timesum_myr = seconds_to_years(timesum_seconds)/1e6

        if timesum_myr >= solver.pelagic_sedimentation_rate_reduction_time
            pelagic_sedimentation_rate /= solver.pelagic_sedimentation_rate_reduction_factor
        end
    end
    print_info("Pelagic sedimentation rate (mm/yr): $(meters_per_seconds_to_mm_per_yr(pelagic_sedimentation_rate))", level=2)
    sediment_thickness_total = copy(solver.sediment_thickness_initial)
    for istep in 1:solver.nsteps
        if solver.use_print_debug
            print_timestep_info(solver, istep)
        end
        topo_grid_diffusivity = make_topo_diffusivity_grid(
            solver.topo_gridx, water_depth_x, downstream_distances_x,
            solver.sediment_transport_parameters,
            solver.use_constant_diffusivity
        )
        (
            topo_grid_pelagic_sedimentation_rate
        ) = calculate_pelagic_sedimentation_rate_grid(
            pelagic_sedimentation_rate,
            solver.topo_gridx, water_depth_x
        )
        solution_array = solve_downhill_diffusion(
            solver.topo_gridx, solver.topo_gridy, topo_grid_diffusivity,
            topo_grid_pelagic_sedimentation_rate,
            solver.basic_grid_x_dimensions, solver.transport_timestep,
            solver.sediment_transport_parameters
        )
        copy_array_1d!(solution_array, solver.topo_gridy)
        if solver.use_compaction_correction
            if solver.compaction_correction_type == "constant_property" || isnothing(model)
                (
                    topo_gridy_corrected, sediment_thickness_total, _
                ) = apply_compaction_correction(
                    solver.sediment_transport_parameters,
                    solver.topo_gridy_initial,
                    solver.topo_gridy,
                    sediment_thickness_total
                )
                copy_array_1d!(topo_gridy_corrected, solver.topo_gridy)
            elseif solver.compaction_correction_type == "variable_property" && !isnothing(model)
                (
                    topo_gridy_corrected, sediment_thickness_total, _
                ) = apply_compaction_correction_for_topography_and_markers(
                        model,
                        solver.sediment_transport_parameters,
                        solver.topo_gridx,
                        solver.topo_gridy_initial,
                        solver.topo_gridy,
                        sediment_thickness_total,
                        markers_indices_sedimentary_basin,
                        markers_indices_sticky,
                        x_sorted_marker_indices_sedimentary_basin,
                        x_sorted_marker_indices_sticky
                    )
                copy_array_1d!(topo_gridy_corrected, solver.topo_gridy)
            end
        end
        copy_array_1d!(solver.topo_gridy, solver.topo_gridy_initial)
        (
            downstream_distances_x, drainage_divides_x
        ) = calculate_downstream_distances_for_nodes(
            solver.topo_gridx, solver.topo_gridy)
        water_depth_x = calculate_water_depth(solver.topo_gridy, sealevel_x)
        update_collections!(
            solver.collections, istep, solver.topo_gridy, 
            drainage_divides_x, water_depth_x
            )
    end

    if solver.use_compaction_correction
        solver.sediment_thickness_total_decompacted = copy(sediment_thickness_total)
        if solver.compaction_correction_type == "constant_property"
            update_compaction_displacement_max!(solver, sediment_thickness_total)
        end
    end
    return nothing
end

function update_compaction_displacement_max!(
    solver::SedimentTransportSolver,
    sediment_thickness_total::Vector{Float64}
)::Nothing
    sediment_thickness_markers_final = calculate_final_thickness_after_burial(
        solver.sediment_transport_parameters.porosity_initial,
        solver.sediment_transport_parameters.depth_decay_term,
        solver.sediment_thickness_initial,
        sediment_thickness_total
    )

    solver.sediment_thickness_initial_compacted = copy(sediment_thickness_markers_final)

    compaction_displacement_max = solver.sediment_thickness_initial - 
        sediment_thickness_markers_final
    solver.compaction_displacement_max = copy(compaction_displacement_max)

    return nothing
end

function calculate_pelagic_sedimentation_rate_grid(
    pelagic_sedimentation_rate::Float64,
    topo_gridx::Vector{Float64},
    water_depth_x::Vector{Float64}
)::Vector{Float64}
    toponum = length(topo_gridx)
    topo_grid_pelagic_sedimentation_rate = zeros(Float64, toponum)
    
    for i in 2:toponum
        water_depth = water_depth_x[i]
        if water_depth > 0.0
            topo_grid_pelagic_sedimentation_rate[i] = pelagic_sedimentation_rate
        else
            topo_grid_pelagic_sedimentation_rate[i] = 0.0
        end
    end
    return topo_grid_pelagic_sedimentation_rate
end

function print_timestep_info(solver::SedimentTransportSolver, istep::Int)::Nothing
    print_info(
        "Working on time step: $istep total transport time (Myr): "
        * "$(seconds_to_years((istep)*solver.transport_timestep)/1e6)", 
        level=2
        )
    return nothing
end

end