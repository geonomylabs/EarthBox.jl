module Interpolation

include("mapping/MarkerGridMapping.jl")
include("weights/WeightFuncs.jl")
include("utils/BilinearAverage.jl")
include("utils/BilinearNumerator.jl")
include("backup/BackupTransportArrays.jl")
include("initialize/InitializeInterpolation.jl")
include("utils/GridToMarker.jl")
include("utils/GridTemperatureToMarkers.jl")
include("utils/heat_funcs/ApplyThermalBC.jl")
include("utils/stokes_funcs/CheckPlasticYield.jl")
include("utils/stokes_funcs/HarmonicAverage.jl")
include("utils/stokes_funcs/FlattenDensity.jl")
include("utils/stokes_funcs/PlasticStrainRateInterpolation.jl")
include("utils/MarkerToHeatArrays.jl")
include("utils/MarkersToStokesArrays.jl")

import EarthBox.ModelDataContainer: ModelData
import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.PrintFuncs: @timeit_memit
import EarthBox.ModelDataContainer: load_parameters!
import EarthBox.InitializationTools: InitializationTools
import EarthBox.UseOptionTools: set_use_option!
import .InitializeInterpolation
import .MarkerToHeatArrays
import .MarkersToStokesArrays
import .HarmonicAverage

const PDATA = get_eb_parameters()

struct ValidInputNames
    iuse_harmonic_avg_normal_viscosity::Symbol
    iuse_initial_temp_interp::Symbol
end

"""
    initialize!(model::ModelData; kwargs...)::Nothing

Initialize interpolation option parameter flags.

# Arguments
- `model::`[`ModelData`](@ref ModelData): The model data container containing the 
    model parameters and arrays.

# Keyword Arguments
- `$(PDATA.iuse_harmonic_avg_normal_viscosity.name)::Int64`:
    - $(PDATA.iuse_harmonic_avg_normal_viscosity.description)
- `$(PDATA.iuse_initial_temp_interp.name)::Int64`:
    - $(PDATA.iuse_initial_temp_interp.description)

"""
function initialize!(
    model::ModelData;
    kwargs...
)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    return nothing
end

function initialize_bilinear_interpolation!(
    model::ModelData, 
    inside_flags::Vector{Int8}
)::Nothing
    @timeit_memit "Finished initializing bilinear interpolation" begin
        # Initialize bilinear interpolation parameters
        # Parameters are calculated that are used for marker-grid interpolation
        # This includes:
        # 1. indices and distance of the upper-left grid node associated with the cell containing the marker
        # 2. marker weights for the four grid nodes and the central node of the cell that contains the marker
        # 3. marker-weight sums for grid nodes that are used in the denominator of the bilinear interpolation calculation
        InitializeInterpolation.initialize_bilinear_interpolation!(model, inside_flags)
    end
    return nothing
end

function interpolate_marker_heat_parameters_to_grid!(
    model::ModelData,
    inside_flags::Vector{Int8}
)::Nothing
    @timeit_memit "Finished interpolating marker heat parameters to grid" begin
        MarkerToHeatArrays.interpolate_markers_to_transport_arrays!(model, inside_flags)
    end
    return nothing
end

function interpolate_marker_stokes_parameters_to_grid!(
    model::ModelData,
    inside_flags::Vector{Int8}
)::Nothing
    @timeit_memit "Finished interpolating marker stokes parameters to grid" begin
        MarkersToStokesArrays.interpolate_markers_to_transport_arrays!(model, inside_flags)
    end
    return nothing
end

function interpolate_basic_viscoplastic_viscosity_to_pressure_grid!(
    model::ModelData
)::Nothing
    # Update the harmonic viscoplastic viscosity array at pressure nodes
    # This array is referred to as etan1 and is calculated using a harmonic
    # average of etas1, the viscoplastic viscosity at basic nodes
    if get_iuse_harmonic(model) == 1
        @timeit_memit "Finished interpolating basic viscoplastic viscosity to pressure grid" begin
            HarmonicAverage.normal_viscosity_from_harmonic_avg_shear_viscosity!(model)
        end
    end
    return nothing
end

function get_iuse_harmonic(model::ModelData)::Int64
    interp_options = model.interpolation.parameters.interp_options
    return interp_options.iuse_harmonic_avg_normal_viscosity.value
end

function interpolate_temperature_back_to_markers_at_startup!(
    model::ModelData,
    inside_flags::Vector{Int8}
)::Nothing
    # Interpolate grid temperature back to markers at startup
    # This method interpolates the grid temperature tk1 back to the marker
    # temperature marker_TK at the start of the simulation if the option is selected
    # This is done to avoid discrepancies between the initial nodal and marker temperatures
    timesum = model.timestep.parameters.main_time_loop.timesum.value
    if get_iuse_initial_temp_interp(model) == 1 && timesum == 0
        @timeit_memit "Finished interpolating temperature back to markers at startup" begin
            GridTemperatureToMarkers.interpolate!(model, inside_flags)
        end
    end
    return nothing
end

function get_iuse_initial_temp_interp(model::ModelData)::Int64
    interp_options = model.interpolation.parameters.interp_options
    return interp_options.iuse_initial_temp_interp.value
end

end # module 