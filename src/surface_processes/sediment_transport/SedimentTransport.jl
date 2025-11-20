module SedimentTransport

include("solver/SedimentTransportSolverManager.jl")

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: load_parameters!
import EarthBox.SedimentThickness: calculate_sediment_thickness_from_markers
import EarthBox.DataStructures: SedimentTransportParameters
import EarthBox.EBCopy: copy_new_topography_to_topography_array
import EarthBox.EBCopy: copy_topography_coordinate_arrays
import EarthBox.Compaction: get_boolean_options_for_compaction
import EarthBox.TimeStep: TimeStepCalculator
import .SedimentTransportSolverManager: SedimentTransportSolver
import .SedimentTransportSolverManager: run_sediment_transport_time_steps!
import .SedimentTransportSolverManager.MarkerAdvection: advect_markers_using_compaction

const PDATA = get_eb_parameters()

struct SedimentTransportModel
    iuse_downhill_diffusion::Union{Int64, Nothing}
    subaerial_slope_diffusivity::Union{Float64, Nothing}
    precipitation_rate::Union{Float64, Nothing}
    subaerial_transport_coefficient::Union{Float64, Nothing}
    submarine_slope_diffusivity::Union{Float64, Nothing}
    submarine_diffusion_decay_depth::Union{Float64, Nothing}
    transport_timestep::Union{Float64, Nothing}
    pelagic_sedimentation_rate::Union{Float64, Nothing}
    pelagic_sedimentation_rate_reduction_factor::Union{Float64, Nothing}
    pelagic_sedimentation_rate_reduction_time::Union{Float64, Nothing}
    iuse_compaction_correction::Union{Int64, Nothing}
end

struct ValidInputNames
    iuse_downhill_diffusion::Symbol
    subaerial_slope_diffusivity::Symbol
    precipitation_rate::Symbol
    subaerial_transport_coefficient::Symbol
    submarine_slope_diffusivity::Symbol
    submarine_diffusion_decay_depth::Symbol
    transport_timestep::Symbol
    pelagic_sedimentation_rate::Symbol
    pelagic_sedimentation_rate_reduction_factor::Symbol
    pelagic_sedimentation_rate_reduction_time::Symbol
    iuse_compaction_correction::Symbol
end

"""
    initialize!(
        model::Union{ModelData, Nothing}; 
        kwargs...
    )::Union{SedimentTransportModel, Nothing}

Initialize sediment transport model parameters.

# Arguments
- `model::`[`ModelData`](@ref ModelData)
    - The model data container containing the model parameters and arrays.

# Keyword Arguments
- `$(PDATA.iuse_downhill_diffusion.name)::Int64`
    - $(PDATA.iuse_downhill_diffusion.description)
- `$(PDATA.subaerial_slope_diffusivity.name)::Float64`
    - $(PDATA.subaerial_slope_diffusivity.description)
- `$(PDATA.precipitation_rate.name)::Float64`
    - $(PDATA.precipitation_rate.description)
- `$(PDATA.subaerial_transport_coefficient.name)::Float64`
    - $(PDATA.subaerial_transport_coefficient.description)
- `$(PDATA.submarine_slope_diffusivity.name)::Float64`
    - $(PDATA.submarine_slope_diffusivity.description)
- `$(PDATA.submarine_diffusion_decay_depth.name)::Float64`
    - $(PDATA.submarine_diffusion_decay_depth.description)
- `$(PDATA.transport_timestep.name)::Float64`
    - $(PDATA.transport_timestep.description)
- `$(PDATA.pelagic_sedimentation_rate.name)::Float64`
    - $(PDATA.pelagic_sedimentation_rate.description)
- `$(PDATA.pelagic_sedimentation_rate_reduction_factor.name)::Float64`
    - $(PDATA.pelagic_sedimentation_rate_reduction_factor.description)
- `$(PDATA.pelagic_sedimentation_rate_reduction_time.name)::Float64`
    - $(PDATA.pelagic_sedimentation_rate_reduction_time.description)
- `$(PDATA.iuse_compaction_correction.name)::Int64`
    - $(PDATA.iuse_compaction_correction.description)
"""
function initialize!(
    model::Union{ModelData, Nothing};
    kwargs...
)::Union{SedimentTransportModel, Nothing}
    if !(model === nothing)
        load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
        return nothing
    else
        return SedimentTransportModel(
            get(kwargs, :iuse_downhill_diffusion, nothing),
            get(kwargs, :subaerial_slope_diffusivity, nothing),
            get(kwargs, :precipitation_rate, nothing),
            get(kwargs, :subaerial_transport_coefficient, nothing),
            get(kwargs, :submarine_slope_diffusivity, nothing),
            get(kwargs, :submarine_diffusion_decay_depth, nothing),
            get(kwargs, :transport_timestep, nothing),
            get(kwargs, :pelagic_sedimentation_rate, nothing),
            get(kwargs, :pelagic_sedimentation_rate_reduction_factor, nothing),
            get(kwargs, :pelagic_sedimentation_rate_reduction_time, nothing),
            get(kwargs, :iuse_compaction_correction, nothing)
        )
    end
end

function run_sediment_transport_model!(
    model::ModelData;
    compaction_correction_type::String="constant_property"
)::Union{Vector{Float64}, Nothing}
    use_compaction_correction = get_boolean_options_for_compaction(model)
    gridx_b = model.grids.arrays.basic.gridx_b.array
    gridt = model.topography.arrays.gridt.array
    topo_gridx, topo_gridy_initial = copy_topography_coordinate_arrays(gridt)
    sediment_transport_parameters = get_sediment_transport_parameters(model)
    y_sealevel = model.topography.parameters.sealevel.y_sealevel.value
    (
        sediment_thickness_markers_initial
    ) = calculate_sediment_thickness_from_markers(model, topo_gridx)
    (
        pelagic_sedimentation_rate
    ) = model.topography.parameters.depo_and_erosion_rates.pelagic_sedimentation_rate.value
    transport_solver = SedimentTransportSolver(
        (gridx_b[1], gridx_b[end]),
        topo_gridx,
        topo_gridy_initial,
        sediment_thickness_markers_initial,
        sediment_transport_parameters,
        y_sealevel,
        pelagic_sedimentation_rate,
        use_collections=false,
        use_print_debug=false,
        use_constant_diffusivity=false,
        use_compaction_correction=use_compaction_correction,
        compaction_correction_type=compaction_correction_type
    )
    run_sediment_transport_time_steps!(transport_solver, model)
    if transport_solver.use_compaction_correction && transport_solver.compaction_correction_type == "constant_property"
        advect_markers_using_compaction(
            model,
            sediment_thickness_markers_initial,
            transport_solver.compaction_displacement_max
        )
    end
    copy_new_topography_to_topography_array(transport_solver.topo_gridy, gridt)
    return transport_solver.sediment_thickness_total_decompacted
end

function get_sediment_transport_parameters(model::ModelData)
    porosity_initial, decay_depth_term = get_sediment_compaction_properties(model)

    downhill_diffusion = model.topography.parameters.downhill_diffusion

    return SedimentTransportParameters(
        downhill_diffusion.subaerial_slope_diffusivity.value,
        downhill_diffusion.precipitation_rate.value,
        downhill_diffusion.subaerial_transport_coefficient.value,
        downhill_diffusion.submarine_slope_diffusivity.value,
        downhill_diffusion.submarine_diffusion_decay_depth.value,
        downhill_diffusion.transport_timestep.value,
        model.timestep.parameters.main_time_loop.timestep.value,
        porosity_initial,
        decay_depth_term
    )
end

function get_sediment_compaction_properties(model::ModelData)
    types = model.materials.dicts.matid_types
    matid_sediment = types["Sediment"][1]
    mat_compaction = model.materials.arrays.mat_compaction.array

    porosity_initial = mat_compaction[matid_sediment, 1]
    decay_depth = mat_compaction[matid_sediment, 2]
    decay_depth_term = 1.0/decay_depth

    return porosity_initial, decay_depth_term
end

end # module 