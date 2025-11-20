module Extraction

include("core/DebugTools.jl")
include("core/Extractable.jl")
include("core/MeltCheck.jl")
include("core/Extracted.jl")
include("core/PartialMeltIndices.jl")
include("core/PartiallyMoltenZone.jl")
include("core/MagmaBody.jl")
include("core/Volcanism.jl")
include("core/MeltVolumetrics.jl")

using Printf
import Plots
import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.ModelDataContainer: ModelData
import EarthBox.PrintFuncs: @timeit_memit, print_info, print_warning, 
    print_melt_extraction_info
import EarthBox.ModelDataContainer: load_parameters!
import EarthBox.ConversionFuncs: celsius_to_kelvin, seconds_to_years
import EarthBox.MagmaFlushState: update_use_extrusion_for_magma_flush
import .MeltVolumetrics
import .Volcanism
import .Extracted
import .Extractable
import .PartialMeltIndices
import .PartiallyMoltenZone
import .MagmaBody
import ..Drainage
import ..MoltenZone
import ..MeltRefraction

const PDATA = get_eb_parameters()

const DEBUG = false

struct ValidInputNames
    iuse_extraction::Symbol
    iuse_gabbroic_fractionation::Symbol
    iuse_shallow_mantle_injection::Symbol
    iuse_random_injection_subdomain::Symbol
    iuse_normal_injection_subdomain::Symbol
    number_of_injection_subdomains::Symbol
    emplacement_temperature::Symbol
    smoothing_radius_drainage::Symbol
    smoothing_radius_fractionation::Symbol
    characteristic_injection_width::Symbol
    mantle_search_width::Symbol
    magma_height_limit::Symbol
    fractionation_threshold_limit::Symbol
    maximum_shallow_injection_depth::Symbol
    extraction_fraction::Symbol
end

"""
    initialize!(model::ModelData; kwargs...)::Nothing

Initialize melt extraction model parameters.

# Arguments
- `model::`[`ModelData`](@ref ModelData): The model data container containing the 
    model parameters and arrays.

# Keyword Arguments

- `$(PDATA.iuse_extraction.name)::Int64`
    - $(PDATA.iuse_extraction.description)
- `$(PDATA.iuse_shallow_mantle_injection.name)::Int64`
    - $(PDATA.iuse_shallow_mantle_injection.description)
- `$(PDATA.smoothing_radius_drainage.name)::Float64`
    - $(PDATA.smoothing_radius_drainage.description)
- `$(PDATA.smoothing_radius_fractionation.name)::Float64`
    - $(PDATA.smoothing_radius_fractionation.description)
- `$(PDATA.characteristic_injection_width.name)::Float64`
    - $(PDATA.characteristic_injection_width.description)
- `$(PDATA.mantle_search_width.name)::Float64`
    - $(PDATA.mantle_search_width.description)
- `$(PDATA.magma_height_limit.name)::Float64`
    - $(PDATA.magma_height_limit.description)
- `$(PDATA.iuse_random_injection_subdomain.name)::Int64`
    - $(PDATA.iuse_random_injection_subdomain.description)
- `$(PDATA.iuse_normal_injection_subdomain.name)::Int64`
    - $(PDATA.iuse_normal_injection_subdomain.description)
- `$(PDATA.number_of_injection_subdomains.name)::Int64`
    - $(PDATA.number_of_injection_subdomains.description)
- `$(PDATA.emplacement_temperature.name)::Float64`
    - $(PDATA.emplacement_temperature.description)
- `$(PDATA.iuse_gabbroic_fractionation.name)::Int64`
    - $(PDATA.iuse_gabbroic_fractionation.description)
- `$(PDATA.fractionation_threshold_limit.name)::Float64`
    - $(PDATA.fractionation_threshold_limit.description)
- `$(PDATA.maximum_shallow_injection_depth.name)::Float64`
    - $(PDATA.maximum_shallow_injection_depth.description)
- `$(PDATA.extraction_fraction.name)::Float64`
    - $(PDATA.extraction_fraction.description)
"""
function initialize!(model::ModelData; kwargs...)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    return nothing
end

""" Update melt extraction model.

This method updates the melt extraction model by calculating melt drainage 
divides, extracting melt from the mantle to magma bodies, calculating 
dimensions of molten domains and transforming markers for melt refraction if 
the extractable melt fraction is less than zero.
"""
function update_melt_extraction!(
    model::ModelData,
    inside_flags::Vector{Int8},
    output_dir::String
)::Nothing
    @timeit_memit "Finished updating melt extraction" begin
        update_melt_drainage_divides!(model, output_dir)
        extract_melt_from_mantle_to_magma_bodies!(model)
        calculate_dimensions_of_molten_domains!(model)
        transform_markers_for_melt_refraction!(model, inside_flags)
    end
    return nothing
end

""" Update melt drainage divides.

This method calculates melt drainage divides and redistribute eruption volumes 
to new drainage basins.
"""
function update_melt_drainage_divides!(
    model::ModelData,
    output_dir::Union{String, Nothing}
)
    iuse_extraction = model.melting.parameters.options.iuse_extraction.value
    iuse_melting = model.melting.parameters.options.iuse_melting.value
    iuse_eruption_interval = model.melting.parameters.extrusion.iuse_eruption_interval.value

    if iuse_melting == 1 && iuse_extraction == 1
        @timeit_memit "Finished calculating melt drainage divides" begin
            topo_gridx, partial_melt_gridy, divides_x = 
                Drainage.calculate_melt_drainage_divides!(model)
        end
        if iuse_eruption_interval == 1
            @timeit_memit "Finished redistributing extrusion volumes to new drainage basins" begin
                Drainage.redistribute_extrusion_volumes_to_new_drainage_basins!(model)
            end
        end
        if DEBUG
            make_debug_plots(model, output_dir, topo_gridx, partial_melt_gridy, divides_x)
        end
    end
end

""" Extract melt from mantle to magma bodies.

This module extracts melt from mantle to magma bodies for all drainage basins.
"""
function extract_melt_from_mantle_to_magma_bodies!(
    model::ModelData
)::Nothing
    iuse_extraction = model.melting.parameters.options.iuse_extraction.value
    iuse_melting = model.melting.parameters.options.iuse_melting.value
    if iuse_melting == 1 && iuse_extraction == 1
        @timeit_memit "Finished extracting melt from mantle to magma bodies" begin
            extract_melt_from_mantle_for_all_drainage_basins!(model)
        end
    end
    return nothing
end

""" Execute melt extraction steps.

This method calculates dimensions of molten domain for all drainage basins.
"""
function calculate_dimensions_of_molten_domains!(
    model::ModelData
)::Nothing
    iuse_extraction = model.melting.parameters.options.iuse_extraction.value
    iuse_melting = model.melting.parameters.options.iuse_melting.value
    
    if iuse_melting == 1 && iuse_extraction == 1
        @timeit_memit "Finished calculating molten zone dimensions" begin
            MoltenZone.call_calculate_dimensions_of_molten_domain_all_drainage!(model)
        end
    end
    return nothing
end

function transform_markers_for_melt_refraction!(
    model::ModelData, 
    inside_flags::Vector{Int8}
)::Nothing
    iuse_extraction = model.melting.parameters.options.iuse_extraction.value
    iuse_melting = model.melting.parameters.options.iuse_melting.value
    if iuse_melting == 1 && iuse_extraction == 1
        @timeit_memit "Finished updating marker type for melt refraction" begin
            MeltRefraction.update_marker_type_for_melt_refraction!(model, inside_flags)
        end
    end
    return nothing
end

function extract_melt_from_mantle_for_all_drainage_basins!(model::ModelData)
    initial_magma_flush_steps = model.melting.parameters.extrusion.initial_magma_flush_steps.value
    magma_flush_factor = model.melting.parameters.extrusion.magma_flush_factor.value
    extract_melt_in_drainage_basins!(model, initial_magma_flush_steps, magma_flush_factor)
end

""" Extract partial melt to magma bodies at the top of partial melt zone.

Each drainage basin is considered separately.

First, find mantle markers with partial melt (i.e. where melt is stable) and 
calculate the total incremental volume of melt to extract. Then calculate the 
total number of mantle particles that need to be converted to magma, extracted 
melt fraction and extractable melt fraction. Finally, for each particle that 
needs to be converted to molten material, loop over all partially molten mantle 
particles to find the shallowest (skip already converted particles) and convert 
this shallowest particle to molten mantle by changing rock id for molten mantle.

The partial_melt_marker_indices array contains a list of particle ID's for 
partially molten mantle material.

This function only extracts melt from the mantle to the magma body. It does not 
include processes occurring in the crust.
"""
function extract_melt_in_drainage_basins!(
    model::ModelData,
    initial_magma_flush_steps::Int,
    magma_flush_factor::Float64
)::Nothing
    iuse_extrusion = model.melting.parameters.extrusion.iuse_extrusion.value
    iuse_eruption_interval = model.melting.parameters.extrusion.iuse_eruption_interval.value
    extraction_fraction = model.melting.parameters.extraction.extraction_fraction.value
    characteristic_injection_width = model.melting.parameters.extraction.characteristic_injection_width.value
    magma_height_limit = model.melting.parameters.extraction.magma_height_limit.value

    avg_shallow_partial_melt_xcoors = model.melting.arrays.extraction.avg_shallow_partial_melt_xcoors.array
    avg_shallow_partial_melt_ycoors = model.melting.arrays.extraction.avg_shallow_partial_melt_ycoors.array

    iuse_extrusion = update_use_extrusion_for_magma_flush(model, iuse_extrusion)

    mantle_melting_mat_ids = get_mantle_melting_ids(model)
    mantle_emplacement_mat_ids = get_mantle_emplacement_ids(model)
    Extractable.update_extractable_meltfrac!(model, mantle_melting_mat_ids)
    Extracted.update_extracted_meltfrac!(model, mantle_melting_mat_ids)

    ndrainage_basin = model.melting.parameters.extraction.ndrainage_basin.value
    xstart_drainage = model.melting.arrays.extraction.xstart_drainage.array
    xend_drainage = model.melting.arrays.extraction.xend_drainage.array

    melt_residuals = model.melting.arrays.extraction.melt_residuals.array
    extrusion_volumes = model.melting.arrays.extraction.extrusion_volumes.array

    ntimestep = model.timestep.parameters.main_time_loop.ntimestep.value
    timesum = model.timestep.parameters.main_time_loop.timesum.value
    timesum_myr = seconds_to_years(timesum)/1e6

    total_melt_volume_check_marker_units = MeltVolumetrics.calculate_melt_volume(
        model, mantle_melting_mat_ids,
        xstart=0.0, xend=model.grids.parameters.geometry.xmax.value
        )

    characteristic_magmatic_crust_height = 
        calculate_characteristic_new_crust_height(model, total_melt_volume_check_marker_units)

    model.melting.parameters.extrusion.characteristic_magmatic_crust_height.value = 
        characteristic_magmatic_crust_height

    print_characteristic_new_crust_height(timesum_myr, characteristic_magmatic_crust_height)

    total_melt_volume_sum_marker_units = 0.0
    for i in 1:ndrainage_basin
        xstart = xstart_drainage[i]
        xend = xend_drainage[i]

        # Define marker indices that are in the partially molten zone associated
        # with the current drainage basin.
        (
            partial_melt_marker_indices, nmarkers_partial_melt
        ) = PartialMeltIndices.get_partial_melt_marker_indices(
                model, xstart=xstart, xend=xend)

        melt_volume_mantle_marker_units = MeltVolumetrics.calculate_melt_volume(
            model, mantle_melting_mat_ids, xstart=xstart, xend=xend)

        total_melt_volume_sum_marker_units += melt_volume_mantle_marker_units

        melt_volume_mantle_marker_units += melt_residuals[i]
        melt_volume_mantle_marker_units *= extraction_fraction

        # nmarkers_magma is the number of markers in the partially molten zone
        # that will be converted into molten magma
        nmarkers_magma_mantle = floor(Int, melt_volume_mantle_marker_units)

        (
            magma_production_rate_m3_yr, avg_marker_volume_m3
        ) = MeltVolumetrics.calculate_magma_production_rate(model, nmarkers_magma_mantle)

        # Calculate the fractional residual incremental melt volume that will
        # be used to determine the number of markers that will be converted to
        # magma during the next time step.
        melt_volume_residual_marker_units = 
            MeltVolumetrics.calculate_residual_melt_volume(
                melt_volume_mantle_marker_units, nmarkers_magma_mantle)

        nmarkers_volcanics = Volcanism.calculate_number_of_volcanic_markers(
            model,
            nmarkers_magma_mantle,
            characteristic_magmatic_crust_height,
            mantle_melting_mat_ids,
            xstart=xstart,
            xend=xend,
            initial_magma_flush_steps=initial_magma_flush_steps,
            magma_flush_factor=magma_flush_factor
            )

        nmarkers_magma_mantle = Volcanism.update_number_of_magma_markers_for_extrusion(
            nmarkers_magma_mantle, nmarkers_volcanics, iuse_extrusion)

        injection_width = calculate_injection_width(
            model,
            characteristic_injection_width,
            magma_height_limit,
            nmarkers_magma_mantle
            )

        (
            xshallow_partial_melt_avg, yshallow_partial_melt_avg
        ) = MagmaBody.extract_partial_melt_and_make_magma_body(
                model,
                mantle_emplacement_mat_ids,
                nmarkers_magma_mantle,
                nmarkers_partial_melt,
                partial_melt_marker_indices,
                injection_width
            )

        avg_shallow_partial_melt_xcoors[i] = xshallow_partial_melt_avg
        avg_shallow_partial_melt_ycoors[i] = yshallow_partial_melt_avg

        # The following block of code was added to avoid problems with magma
        # flushing if serpentinization is used. How serpentinization causes
        # magma flushing is not clear.
        if Volcanism.magma_flush_event(ntimestep, initial_magma_flush_steps)
            extrusion_volume_m3 = 0.0
        else
            extrusion_volume_m3 = Volcanism.calculate_extrusion_volume(
                model, nmarkers_volcanics, iuse_extrusion)
        end

        melt_residuals[i] = melt_volume_residual_marker_units

        # If an eruption interval is used material available for extrusion is 
        # summed over multiple time steps.
        if iuse_eruption_interval == 0
            extrusion_volumes[i] = extrusion_volume_m3
        else
            extrusion_volumes[i] += extrusion_volume_m3
        end

        print_info = true
        if print_info
            print_extraction_info(
                i, timesum_myr, nmarkers_partial_melt, 
                melt_volume_mantle_marker_units,
                nmarkers_magma_mantle, melt_volume_residual_marker_units,
                nmarkers_volcanics, extrusion_volume_m3,
                xstart, xend, characteristic_injection_width, injection_width,
                xshallow_partial_melt_avg,
                magma_production_rate_m3_yr, avg_marker_volume_m3)
        end
    end

    print_melt_volume_check(
        model, total_melt_volume_sum_marker_units,
        total_melt_volume_check_marker_units
    )
    return nothing
end

""" Calculate characteristic new crust height.

The characteristic new crust height is the height of new crust that would be 
formed if all the melt in the mantle was extracted and emplaced at the top of 
the model in a column with a width equal to the full extension velocity times 
the time step.
"""
function calculate_characteristic_new_crust_height(
    model::ModelData,
    total_incremental_mantle_melt_volume_marker_units::Float64
)::Float64
    total_incremental_mantle_melt_volume_m3 = 
        total_incremental_mantle_melt_volume_marker_units *
        MeltVolumetrics.calculate_avg_marker_volume(model)
    
    full_extension_velocity_m_s = model.bcs.parameters.velocity.full_velocity_extension.value
    timestep_seconds = model.timestep.parameters.main_time_loop.timestep.value
    characteristic_new_crust_width_m = full_extension_velocity_m_s * timestep_seconds
    characteristic_new_crust_area_m2 = characteristic_new_crust_width_m * 1.0

    if characteristic_new_crust_area_m2 > 0.0
        characteristic_new_crust_height_m = total_incremental_mantle_melt_volume_m3 / characteristic_new_crust_area_m2
    else
        characteristic_new_crust_height_m = 0.0
    end
    
    return characteristic_new_crust_height_m
end

function print_melt_volume_check(
    model::ModelData,
    total_melt_volume_sum_marker_units::Float64,
    total_melt_volume_check_marker_units::Float64
)
    avg_marker_volume_m3 = MeltVolumetrics.calculate_avg_marker_volume(model)
    msg = @sprintf(
        "total melt volume (all basins) (m^3): %.2f : melt volume from mantle (m^3) : %.2f",
        total_melt_volume_sum_marker_units * avg_marker_volume_m3,
        total_melt_volume_check_marker_units * avg_marker_volume_m3
    )
    print_info(msg, level=2)
end

""" Calculate injection width.

# Arguments
- model::ModelData
    - Model data structure.
- characteristic_injection_width::Float64
    - Characteristic injection width in meters.
- magma_height_limit::Float64
    - Maximum magma height in meters.
- nmarkers_magma::Float64
    - Number of magma markers and the volume of magma in marker units.

# Returns
- injection_width::Float64
    - Injection width in meters.
"""
function calculate_injection_width(
    model::ModelData,
    characteristic_injection_width::Float64,
    magma_height_limit::Float64,
    nmarkers_magma::Int64
)::Float64
    mxstep = model.markers.parameters.distribution.mxstep.value
    mystep = model.markers.parameters.distribution.mystep.value
    volume_per_marker_m3 = mxstep * mystep
    magma_volume_m3 = nmarkers_magma * volume_per_marker_m3
    magma_height = magma_volume_m3 / characteristic_injection_width
    
    if magma_height > magma_height_limit
        injection_width = magma_volume_m3 / magma_height_limit
    else
        injection_width = characteristic_injection_width
    end
    
    return injection_width
end

function print_characteristic_new_crust_height(
    timesum_myr::Float64,
    characteristic_new_crust_height::Float64
)
    msg = @sprintf(
        "Characteristic new crust height (meters): %.2f model time (Myr): %.2f",
        characteristic_new_crust_height, timesum_myr
    )
    print_info(msg, level=2)
end

function print_extraction_info(
    i::Int,
    timesum_myr::Float64,
    nmarkers_partial_melt::Int,
    melt_volume_mantle_marker_units::Float64,
    nmarkers_magma_mantle::Int,
    melt_volume_residual_marker_units::Float64,
    nmarkers_volcanics::Int,
    extrusion_volume_m3::Float64,
    xstart::Float64,
    xend::Float64,
    characteristic_injection_width::Float64,
    injection_width::Float64,
    xshallow_partial_melt_avg::Float64,
    magma_production_rate_m3_yr::Float64,
    avg_marker_volume_m3::Float64
)
    msg = @sprintf(
        "Melt extraction info for drainage basin: %d", i
    )
    print_melt_extraction_info(msg, level=2)
    msg = @sprintf(
        "xstart, xend: %.2f %.2f timesum (Myr): %.2f", 
        xstart, xend, timesum_myr
    )
    print_melt_extraction_info(msg, level=3)
    msg = @sprintf(
        "nmarkers_partial_melt: %d timesum (Myr): %.2f", 
        nmarkers_partial_melt, timesum_myr
    )
    print_melt_extraction_info(msg, level=3)
    msg = @sprintf(
        "extractable melt volume in mantle (m^3): %.2f timesum (Myr): %.2f",
        melt_volume_mantle_marker_units * avg_marker_volume_m3, timesum_myr
    )
    print_melt_extraction_info(msg, level=3)
    msg = @sprintf(
        "melt volume in mantle for intrusion (volcanics removed) (m^3): %.2f timesum (Myr): %.2f",
        nmarkers_magma_mantle * avg_marker_volume_m3, timesum_myr
    )
    print_melt_extraction_info(msg, level=3)
    msg = @sprintf(
        "residual melt volume in mantle (m^3): %.2f timesum (Myr): %.2f",
            melt_volume_residual_marker_units * avg_marker_volume_m3, timesum_myr)
    print_melt_extraction_info(msg, level=3)
    msg = @sprintf(
        "number of volcanic markers: %d timesum (Myr): %.2f",
        nmarkers_volcanics, timesum_myr
    )
    print_melt_extraction_info(msg, level=3)
    msg = @sprintf(
        "extrusion volume (m^3): %.2f timesum (Myr): %.2f",
            extrusion_volume_m3, timesum_myr)
    print_melt_extraction_info(msg, level=3)
    msg = @sprintf(
        "magma characteristic injection width (m): %.2f timesum (Myr): %.2f",
            characteristic_injection_width, timesum_myr)
    print_melt_extraction_info(msg, level=3)
    msg = @sprintf(
        "magma injection width (m): %.2f timesum (Myr): %.2f",
            injection_width, timesum_myr)
    print_melt_extraction_info(msg, level=3)
    msg = @sprintf(
        "average x-coor of shallowest partially molten mantle marker (m): %.2f timesum (Myr): %.2f",
            xshallow_partial_melt_avg, timesum_myr)
    print_melt_extraction_info(msg, level=3)
    msg = @sprintf(
        "magma production rate (m^3/yr): %.2f timesum (Myr): %.2f",
            magma_production_rate_m3_yr, timesum_myr)
    print_melt_extraction_info(msg, level=3)
end

"""
    get_mantle_melting_ids(model::ModelData)::Vector{Int16}

Get all material IDs of mantle melting materials.

Returns
-------
mantle_melting_mat_ids: Vector{Int16}
    Material IDs in the mantle that can undergo melting.
"""
function get_mantle_melting_ids(model::ModelData)::Vector{Int16}
    matid_types = model.materials.dicts.matid_types

    fertile = matid_types["UltramaficMantleFertile"]
    nfertile = length(fertile)

    partially_molten = matid_types["UltramaficMantlePartiallyMolten"]
    npartially_molten = length(partially_molten)

    refractory = matid_types["UltramaficMantleRefactory"]
    nrefractory = length(refractory)

    ntotal = nfertile + npartially_molten + nrefractory

    mantle_melting_mat_ids = zeros(Int16, ntotal)
    icount = 1
    for i in 1:nfertile
        mantle_melting_mat_ids[icount] = fertile[i]
        icount += 1
    end
    for i in 1:npartially_molten
        mantle_melting_mat_ids[icount] = partially_molten[i]
        icount += 1
    end
    for i in 1:nrefractory
        mantle_melting_mat_ids[icount] = refractory[i]
        icount += 1
    end
    return mantle_melting_mat_ids
end

""" Get all material IDs of mantle emplacement materials.

This was added to allow melt emplacement within serpentinite.

# Arguments
- model::ModelData
    - Model data structure.

- use_serpentinite::Bool
    - Boolean flag that activates the use of serpentinite.

# Returns
- mantle_emplacement_mat_ids::Vector{Int16}
    - Material IDs in the mantle where melt can be emplaced.
"""
function get_mantle_emplacement_ids(
    model::ModelData,
    use_serpentinite::Bool=false
)::Vector{Int16}
    matid_types = model.materials.dicts.matid_types

    fertile = matid_types["UltramaficMantleFertile"]
    nfertile = length(fertile)

    partially_molten = matid_types["UltramaficMantlePartiallyMolten"]
    npartially_molten = length(partially_molten)

    refractory = matid_types["UltramaficMantleRefactory"]
    nrefractory = length(refractory)

    serpentinite = matid_types["Serpentinite"]
    nserpentinite = length(serpentinite)

    ntotal = nfertile + npartially_molten + nrefractory

    if use_serpentinite
        ntotal += nserpentinite
    end

    mantle_emplacement_mat_ids = zeros(Int16, ntotal)
    icount = 1
    for i in 1:nfertile
        mantle_emplacement_mat_ids[icount] = fertile[i]
        icount += 1
    end
    for i in 1:npartially_molten
        mantle_emplacement_mat_ids[icount] = partially_molten[i]
        icount += 1
    end
    for i in 1:nrefractory
        mantle_emplacement_mat_ids[icount] = refractory[i]
        icount += 1
    end
    if use_serpentinite
        for i in 1:nserpentinite
            mantle_emplacement_mat_ids[icount] = serpentinite[i]
            icount += 1
        end
    end
    return mantle_emplacement_mat_ids
end

""" Make plot of topography and divides.
"""
function make_debug_plots(
    model::ModelData,
    output_dir::Union{String, Nothing},
    topo_gridx::Vector{Float64},
    partial_melt_gridy::Vector{Float64},
    divides_x::Vector{Float64}
)
    println("")
    println(">> Making plot of melt drainage divides")
    ntimestep = model.timestep.parameters.main_time_loop.ntimestep.value
    ymin = 0.0
    ymax = 160_000.0
    
    # print min and max of partial melt gridy
    println(">> Min partial melt gridy: ", minimum(partial_melt_gridy))
    println(">> Max partial melt gridy: ", maximum(partial_melt_gridy))
    
    dpi = 150
    figsize = (5, 5)
    figsize_pixels = (figsize[1] * dpi, figsize[2] * dpi)

    p = Plots.plot(
        topo_gridx, partial_melt_gridy, label="Top of Partial Melt Zone", 
        color=:blue, linewidth=2.0, linestyle=:solid,
        aspect_ratio=:auto, dpi=dpi, size=figsize_pixels
        )
    Plots.plot!(
        p, divides_x, ones(length(divides_x)) * ymin,
        seriestype=:scatter, label="Drainage Divides", color=:green,
        marker=:circle, markersize=4.0, markerstrokecolor=:transparent, 
        markerstrokewidth=0.0
        )

    Plots.xlabel!("x (m)")
    Plots.ylabel!("Partial Melt Topography")
    Plots.title!("Melt drainage divides")
    Plots.ylims!(ymin, ymax)  # Set y-axis limits
    Plots.yaxis!(:flip)
    
    plot_name = "melt_drainage_divides" * string(ntimestep) * ".png"
    filepath = joinpath(output_dir, plot_name)
    println(">> Saving drainage plot to: ", filepath)
    Plots.savefig(p, filepath)
end

function make_derivative_plot!(
    model::ModelData,
    output_dir::Union{String, Nothing},
    topo_gridx::Vector{Float64},
    partial_melt_gridy::Vector{Float64},
    divides_x::Vector{Float64}
)
    println("")
    println(">> Making plot of melt drainage divides")
    ymin = -0.5
    ymax = 0.5
    ntimestep = model.timestep.parameters.main_time_loop.ntimestep.value
    derivative = gradient(partial_melt_gridy, topo_gridx)
    
    # print min and max derivatives
    println(">> Min derivative: ", minimum(derivative))
    println(">> Max derivative: ", maximum(derivative))
    
    # print derivatives around mid point
    mid = div(length(derivative), 2)
    println(">> Derivative at mid point: ", derivative[mid])
    println(">> Derivative at mid point - 1: ", derivative[mid-1])
    println(">> Derivative at mid point + 1: ", derivative[mid+1])
    
    p = Plots.plot(
        topo_gridx, derivative,
        label="Derivative of Top of Partial Melt Zone", color=:blue
        )
    Plots.plot!(
        p, divides_x, ones(length(divides_x)) * ymin,
        seriestype=:scatter, label="Drainage Divides", color=:green,
        marker=:circle, markersize=4.0, markerstrokecolor=:transparent, 
        markerstrokewidth=0.0
    )
    Plots.xlabel!("x (m)")
    Plots.ylabel!("Derivative of Partial Melt Topography")
    Plots.title!("Melt drainage divides")
    Plots.legend()
    Plots.ylims!(ymin, ymax)  # Set y-axis limits
    Plots.yaxis!(:flip)
    
    plot_name = "melt_drainage_derivatives" * string(ntimestep) * ".png"
    filepath = joinpath(output_dir, plot_name)
    println(">> Saving drainage plot to: ", filepath)
    Plots.savefig(p, filepath)
end 

end # module