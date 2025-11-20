module Fractionation

using Plots
import Printf: @sprintf
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ConversionFuncs: seconds_to_years
import EarthBox.SurfaceProcesses: calculate_age_ma
import EarthBox.MathTools: linear_interp_at_x_location
import EarthBox.ModelStructureManager.TopAndBottom: calculate_top_and_bottom_of_layer_opt
import EarthBox.ModelStructureManager.TopAndBottom: calculate_search_radius
import EarthBox.ModelStructureManager.SmoothSurface: smooth_surface
import ..Extraction.MagmaBody: transform_marker_to_magma
import ..Drainage: calculate_top_of_mantle_partial_melt_domain

""" Call loop function to make fractionated gabbroic magma.

This function calls a loop function loops over all markers and transforms 
gabbroic magma marker to layered gabbroic magma if they are within a certain 
distance from the oceanic Moho by updating marker array `marker_matid` and 
resetting several other marker arrays due to the change in composition.

# Arguments
- `model::ModelData`: The model data container
- `output_dir::String`: The output directory for saving debugging plots
- `debug::Bool`: Boolean flag that activates debugging mode. Debugging mode will
    generate plots of oceanic Moho used for fractionation, the smoothed 
    top of the partial melt domain in the mantle and topography.
"""
function make_fractionated_gabbroic_magma!(
    model::ModelData,
    output_dir::String,
    debug::Bool=false
)::Nothing
    fractionation_threshold_distance = 
        model.melting.parameters.extraction.fractionation_threshold_limit.value
    matid_types = model.materials.dicts.matid_types
    matid_gabbroic_magma = matid_types["ExtractedGabbroicMagma"][1]
    matid_layered_gabbroic_magma = matid_types["ExtractedLayeredGabbroicMagma"][1]

    if matid_gabbroic_magma != -1 && matid_layered_gabbroic_magma != -1
        topo_gridx, topo_gridy, oceanic_moho_gridy, partial_melt_gridy = 
            make_fractionated_gabbroic_magma_loop(
                model, matid_gabbroic_magma, matid_layered_gabbroic_magma,
                fractionation_threshold_distance
            )
        if debug
            plot_topo_moho_and_partial_melting_surface(
                model, topo_gridx, topo_gridy, oceanic_moho_gridy,
                partial_melt_gridy, output_dir
            )
        end
    end
    return nothing
end

""" Make fractionated gabbroic magma with fractionated gabbro composition.

This function loops over all markers and transforms gabbroic magma markers
to layered gabbroic magma if they are within a certain distance from the
Moho. A function is called within the loop that updates marker array 
`marker_matid` and resets several other marker arrays due to the change in
composition.

# Arguments
- `model::ModelData`: The model data container
- `matid_gabbroic_magma::Int`: Material id of the gabbroic magma
- `matid_layered_gabbroic_magma::Int`: Material id of the layered gabbroic magma
- `fractionation_threshold_distance::Float64`: The distance from the oceanic Moho 
    within which gabbroic magma markers are transformed to layered gabbroic magma

# Returns
- `topo_gridx`: The topography grid x array
- `topo_gridy`: The topography grid y array
- `moho_gridy`: The oceanic Moho grid y array
- `partial_melt_gridy`: The partial melt grid y array
"""
function make_fractionated_gabbroic_magma_loop(
    model::ModelData,
    matid_gabbroic_magma::Int16,
    matid_layered_gabbroic_magma::Int16,
    fractionation_threshold_distance::Float64
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}
    gridt = model.topography.arrays.gridt.array
    topo_gridx = copy(gridt[1, :])
    topo_gridy = copy(gridt[2, :])

    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marker_matid = model.markers.arrays.material.marker_matid.array

    age_ma = calculate_age_ma(model)

    nsmooth = calculate_nsmooth(model)
    moho_gridy = calculate_oceanic_moho(model, nsmooth=nsmooth)

    use_partial_melt_limit = true
    if use_partial_melt_limit
        _, _, partial_melt_gridy = calculate_top_of_mantle_partial_melt_domain(model)
        # Make sure that the Moho is not below the partial melt domain.
        # This is necessary to avoid the case where gabbroic particles are
        # located below the Moho due to oceanic crustal thinning.
        moho_gridy = min.(moho_gridy, partial_melt_gridy)
    end

    nmarkers = length(marker_x)
    Threads.@threads for imarker in 1:nmarkers
        @inbounds matid = marker_matid[imarker]
        if matid == matid_gabbroic_magma
            @inbounds begin
                x_marker = marker_x[imarker]
                y_marker = marker_y[imarker]
            end
            y_moho = linear_interp_at_x_location(x_marker, topo_gridx, moho_gridy)
            dist_to_moho = abs(y_marker - y_moho)
            if dist_to_moho < fractionation_threshold_distance
                transform_marker_to_magma(model, imarker, age_ma, matid_layered_gabbroic_magma)
            end
        end
    end

    return topo_gridx, topo_gridy, moho_gridy, partial_melt_gridy
end

function calculate_nsmooth(model::ModelData)
    smoothing_radius = 
        model.melting.parameters.extraction.smoothing_radius_fractionation.value
    gridt = model.topography.arrays.gridt.array
    dx_grid = gridt[1, 2] - gridt[1, 1]
    nsmooth = floor(Int, smoothing_radius/dx_grid)
    nsmooth = max(2, nsmooth)
    return nsmooth
end

function calculate_oceanic_moho(model::ModelData; nsmooth::Int=20)
    gridt = model.topography.arrays.gridt.array
    topo_gridx = copy(gridt[1, :])
    topo_gridy = copy(gridt[2, :])
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marker_matid = model.markers.arrays.material.marker_matid.array
    matids_oc = get_matids_oceanic_crust(model)
    mxstep = model.markers.parameters.distribution.mxstep.value
    marker_search_factor = 
        model.topography.parameters.topo_grid.marker_search_factor.value
    search_radius = calculate_search_radius(mxstep, topo_gridx, marker_search_factor)
    _top_oc, bottom_oc = calculate_top_and_bottom_of_layer_opt(
        matids_oc, marker_matid, marker_x,
        marker_y, topo_gridx, search_radius, use_smoothing=false
    )
    set_zero_values_to_topoy(bottom_oc, topo_gridy)
    oceanic_moho_gridy = smooth_surface(bottom_oc, nsmooth=nsmooth)
    return oceanic_moho_gridy
end

""" Calculate Moho depth.

# Arguments
- `model::ModelData`: The model data container
- `nsmooth::Int`: The number of topography nodes over which to calculate a running
    average
"""
function calculate_moho(model::ModelData, nsmooth::Int=20)
    gridt = model.topography.arrays.gridt.array
    topo_gridx = copy(gridt[1, :])
    topo_gridy = copy(gridt[2, :])
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marker_matid = model.markers.arrays.material.marker_matid.array
    matids_mantle = get_matids_mantle(model)
    mxstep = model.markers.parameters.distribution.mxstep.value
    marker_search_factor = 
        model.topography.parameters.topo_grid.marker_search_factor.value
    search_radius = calculate_search_radius(mxstep, topo_gridx, marker_search_factor)
    top_mantle, _bottom_mantle = calculate_top_and_bottom_of_layer_opt(
        matids_mantle, marker_matid, marker_x,
        marker_y, topo_gridx, search_radius, use_smoothing=false
    )
    set_zero_values_to_topoy(top_mantle, topo_gridy)
    moho_gridy = smooth_surface(top_mantle, nsmooth=nsmooth)
    return moho_gridy
end

function get_matids_mantle(model::ModelData)
    matid_types = model.materials.dicts.matid_types
    ids_mantle = vcat(
        matid_types["UltramaficMantleFertile"],
        matid_types["UltramaficMantlePartiallyMolten"],
        matid_types["UltramaficMantleRefactory"]
    )
    return ids_mantle
end

function get_matids_oceanic_crust(model::ModelData)
    matid_types = model.materials.dicts.matid_types
    ids_oceanic_crust = vcat(
        matid_types["SolidifiedGabbro"],
        matid_types["SolidifiedGabbroPartiallyMolten"],
        matid_types["SolidifiedLayeredGabbro"],
        matid_types["SolidifiedLayeredGabbroPartiallyMolten"],
        matid_types["ExtractedLayeredGabbroicMagma"],
        matid_types["ExtrudedGabbroicMagma"],
        matid_types["SolidifiedBasalt"],
        matid_types["Sediment"]
    )
    return ids_oceanic_crust
end

""" Set zero values to topoy.

# Arguments
- `array`: The array to be updated with ysize at zero values
- `topo_gridy`: The topography grid y array
"""
function set_zero_values_to_topoy(array::Vector{Float64}, topo_gridy::Vector{Float64})
    xnum = length(array)
    Threads.@threads for i in 1:xnum
        @inbounds begin
            if array[i] == 0.0
                array[i] = topo_gridy[i]
            end
        end
    end
end

function plot_topo_moho_and_partial_melting_surface(
    model::ModelData,
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64},
    moho_gridy::Vector{Float64},
    partial_melt_gridy::Vector{Float64},
    output_dir::String
)::Nothing
    ntimestep = model.timestep.parameters.main_time_loop.ntimestep.value
    timesum = model.timestep.parameters.main_time_loop.timesum.value
    timesum_yr = seconds_to_years(timesum)
    timesum_myr = timesum_yr/1e6
    timesum_str = @sprintf("%.2fMyr", timesum_myr)
    
    dpi = 150
    figsize = (5, 5)
    figsize_pixels = figsize .* dpi
    p = plot(
        topo_gridx, topo_gridy, label="Topography", 
        color=:black, linewidth=2, linestyle=:solid,
        size=figsize_pixels
    )
    plot!(
        p, topo_gridx, moho_gridy, label="Moho", color=:red,
        linewidth=2, linestyle=:solid
        )
    plot!(
        p, topo_gridx, partial_melt_gridy, label="Partial melt", color=:blue,
        linewidth=2, linestyle=:solid
    )
    xlabel!("Distance (m)")
    ylabel!("Depth (m)")
    ylims!(0, 50_000)
    yflip!(true)
    
    filepath = joinpath(output_dir, "crust_$(ntimestep)_$(timesum_str).png")
    println(">> Saving crust plot to $(filepath)")
    savefig(p, filepath)
    return nothing
end

end # module 