module MeltModel

include("utils/GetExtractedMoltenIDs.jl")
include("drainage/Drainage.jl")
include("molten_zone/MoltenZone.jl")
include("melt_refraction/MeltRefraction.jl")
include("extraction/Extraction.jl")
include("extrusion/Extrusion.jl")
include("fractionation/Fractionation.jl")
include("melt_fraction/MeltFraction.jl")
include("melt_damage/MeltDamage.jl")
include("melt_properties/MeltPropertiesOpt.jl")
include("melt_rheology/MeltRheology.jl")
include("partial_melting/PartialMelting.jl")
include("solidification/Solidification.jl")
include("rock_props/MeltRockProps.jl")
include("utils/Debug.jl")

using Profile
import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.PrintFuncs: @timeit_memit
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: load_parameters!
import EarthBox.ConversionFuncs: kelvin_to_celsius
import .MeltFraction
import .PartialMelting
import .MeltRheology
import .Extrusion
import .Extraction
import .Solidification
import .Fractionation
import .MeltDamage

const PDATA = get_eb_parameters()

const DEBUG = false

struct ValidInputNames
    iuse_melting::Symbol
    iuse_melt_viscosity::Symbol
    iuse_melt_thermal_props::Symbol
    iuse_depletion_density::Symbol
    iuse_exponential_viscosity_reduction::Symbol
    viscosity_melt::Symbol
end

"""
    initialize!(model::ModelData; kwargs...)::Nothing

Initialize melt model parameters.

# Arguments
- `model::`[`ModelData`](@ref ModelData): The model data container containing the 
    model parameters and arrays.

# Keyword Arguments

- `$(PDATA.iuse_melting.name)::Int64`
    - $(PDATA.iuse_melting.description)
- `$(PDATA.iuse_melt_viscosity.name)::Int64`
    - $(PDATA.iuse_melt_viscosity.description)
- `$(PDATA.viscosity_melt.name)::Float64`
    - $(PDATA.viscosity_melt.description)
- `$(PDATA.iuse_exponential_viscosity_reduction.name)::Int64`
    - $(PDATA.iuse_exponential_viscosity_reduction.description)
- `$(PDATA.iuse_melt_thermal_props.name)::Int64`
    - $(PDATA.iuse_melt_thermal_props.description)
- `$(PDATA.iuse_depletion_density.name)::Int64`
    - $(PDATA.iuse_depletion_density.description)

"""
function initialize!(model::ModelData; kwargs...)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    return nothing
end

"""
    update_marker_melt_model(model::MeltModel, pymodel::ModelData, output_dir::String)

Update marker melt model for the current timestep by transforming markers for
gabbroic fractionation, updating the marker melt fraction array and transforming
markers for solidification and partial melting processes.
"""
function update_marker_melt_model!(
    model::ModelData,
    inside_flags::Vector{Int8},
    output_dir::String
)::Nothing
    transform_markers_for_gabbroic_fractionation!(model, output_dir)
    update_marker_melt_fraction!(model, inside_flags)
    transform_markers_for_melt_processes!(model, inside_flags)
    return nothing
end

"""
Transform markers for gabbroic fractionation by updating the marker material id
array to reflect the fractionated gabbroic melt composition if marker location
is within a threshold distance from the oceanic Moho.
"""
function transform_markers_for_gabbroic_fractionation!(
    model::ModelData,
    output_dir::String
)::Nothing
    timesum = model.timestep.parameters.main_time_loop.timesum.value
    iuse_melting = model.melting.parameters.options.iuse_melting.value
    if timesum > 0 && iuse_melting == 1
        iuse_gabbroic_fractionation = 
            model.melting.parameters.options.iuse_gabbroic_fractionation.value
        if iuse_gabbroic_fractionation == 1
            @timeit_memit "Finished transforming markers for gabbroic fractionation" begin
                Fractionation.make_fractionated_gabbroic_magma!(model, output_dir, DEBUG)
            end
        end
    end 
    return nothing
end

function update_marker_melt_fraction!(
    model::ModelData,
    inside_flags::Vector{Int8}
)::Nothing
    timesum = model.timestep.parameters.main_time_loop.timesum.value
    iuse_melting = model.melting.parameters.options.iuse_melting.value
    if timesum > 0 && iuse_melting == 1
        @timeit_memit "Finished updating marker melt fraction" begin
            MeltFraction.update_marker_melt_fraction!(model, inside_flags)
        end
    end
    return nothing
end

"""
Transform marker composition based on current temperature and pressure by updating
the marker material id array to reflect solidification and partial melting.
"""
function transform_markers_for_melt_processes!(
    model::ModelData,
    inside_flags::Vector{Int8}
)::Nothing
    transform_markers_for_solidification!(model)
    transform_markers_for_partial_melting!(model, inside_flags)
    return nothing
end

"""
Transform marker composition of purely molten marker composition for
solidification if temperature falls below the liquidus.
"""
function transform_markers_for_solidification!(
    model::ModelData
)::Nothing
    timesum = model.timestep.parameters.main_time_loop.timesum.value
    iuse_melting = model.melting.parameters.options.iuse_melting.value
    if timesum > 0 && iuse_melting == 1
        @timeit_memit "Finished applying melt solidification model" begin
            Solidification.solidify!(model)
        end
    end
    return nothing
end

"""
Transform marker composition for partial melting if temperature is between the
solidus and liquidus.
"""
function transform_markers_for_partial_melting!(
    model::ModelData,
    inside_flags::Vector{Int8}
)::Nothing
    timesum = model.timestep.parameters.main_time_loop.timesum.value
    iuse_melting = model.melting.parameters.options.iuse_melting.value
    if timesum > 0 && iuse_melting == 1
        @timeit_memit "Finished updating marker type for partial melting" begin
            PartialMelting.update_solidified_marker_type_for_partial_melting!(model, inside_flags)
        end
    end
    return nothing
end

"""
Update marker rheology for melting by updating marker viscoplastic viscosity and
if partial melt is present reset marker total strain, marker plastic strain and
marker strain rate ratio for melt content.
"""
function update_marker_viscoplastic_viscosity_for_melting!(
    model::ModelData,
    inside_flags::Vector{Int8}
)::Nothing
    timesum = model.timestep.parameters.main_time_loop.timesum.value
    iuse_melting = model.melting.parameters.options.iuse_melting.value
    iuse_melt_viscosity = model.melting.parameters.options.iuse_melt_viscosity.value
    if timesum > 0 && iuse_melting == 1 && iuse_melt_viscosity == 1
        @timeit_memit "Finished updating marker viscoplastic viscosity for melting" begin
            MeltRheology.update_marker_melt_rheology!(model, inside_flags)
        end
    end
    return nothing
end

function update_marker_melt_damage!(model::ModelData)::Nothing
    iuse_melting = model.melting.parameters.options.iuse_melting.value
    if iuse_melting == 1
        MeltDamage.update_melt_damage!(model)
    end
    return nothing
end

""" Update rock properties for melt fraction and depletion for all markers.

Update marker rock properties for melt fraction and depletion by
updating arrays ``marker_rho`` and ``marker_rhocp``. Include the effects
of latent heat in the heat capacity and adiabatic coefficient terms by 
updating marker arrays ``marker_rhocp`` and ``marker_ha``, respectively.

"""
function update_rock_properties_for_melt!(model::ModelData, inside_flags::Vector{Int8})::Nothing
    if use_melt_property_model(model) === true
        @timeit_memit "Finished updating marker rock properties for melt" begin
            MeltRockProps.update_for_melt!(model, inside_flags)
        end
    end
    return nothing
end

function use_melt_property_model(model::ModelData)::Bool
    timesum = model.timestep.parameters.main_time_loop.timesum.value
    iuse_melting = model.melting.parameters.options.iuse_melting.value
    return iuse_melting == 1 && timesum > 0
end

function print_properties_for_max_deltas(
    marker_props::Array{Float64,2}
)::Nothing
    print_delta_info(
        marker_props, 2,
        "Maximum effective heat capacity change due to latent heating/cooling",
        "J/kg/K"
    )
    print_delta_info(
        marker_props, 3,
        "Maximum effective expansivity change due to latent heating/cooling",
        "1/K"
    )
    return nothing
end

function print_delta_info(
    marker_props::Array{Float64,2},
    type_index::Int,
    max_description::String,
    units::String
)::Nothing
    max_value = maximum(marker_props[:,type_index])
    max_index = argmax(marker_props[:,type_index])
    temperature_celsius = kelvin_to_celsius(marker_props[max_index,1])
    pressure_gpa = marker_props[max_index,2]/1e9
    println(
        ">> ", max_description, " ",
        "due to melt fraction: ", max_value, " ", units, " at ",
        "temperature ", temperature_celsius, " C and pressure ", pressure_gpa, " Pa"
    )
    return nothing
end

end # module 