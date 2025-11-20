module SurfaceProcesses

include("utils/CheckMagmaFlush.jl")
include("utils/TopoInterpolation.jl")
include("marker_transformation/MarkerTransformation.jl")
include("sealevel/Sealevel.jl")
include("sediment_transport/SedimentTransport.jl")
include("topography/Topography.jl") # Moved this below sediment transport manager
include("lava_flow/LavaFlowManager.jl")
include("salt/SaltDeposition.jl")

import Profile
import EarthBox.ModelDataContainer: ModelData
import EarthBox.PrintFuncs: @timeit_memit
import EarthBox.Compaction.UpdateProperties: reset_compaction_properties!
import EarthBox.MagmaFlushState: update_use_extrusion_for_magma_flush
import .SedimentTransport: run_sediment_transport_model!
import .LavaFlowManager: lava_flow_manager
import .MarkerTransformation
import .TopoInterpolation
import .SaltDeposition
import .Topography: get_iuse_topo
import .CheckMagmaFlush: check_for_initial_magma_flush


function update_topography!(model::ModelData, inside_flags::Vector{Int8})::Nothing
    Topography.update_topography_using_velocity_field!(model, inside_flags)
    (
        sediment_and_flow_thickness_total_decompacted
    ) = update_topography_using_sediment_transport!(model)
    update_topography_using_lava_flow!(model, sediment_and_flow_thickness_total_decompacted)
    return nothing
end

function update_topography_using_sediment_transport!(
    model::ModelData
)::Union{Vector{Float64}, Nothing}
    iuse_downhill_diffusion = model.obj_dict["iuse_downhill_diffusion"].value
    use_magma_flush = CheckMagmaFlush.check_for_initial_magma_flush(model)
    sediment_and_flow_thickness_total_decompacted = nothing
    if get_iuse_topo(model) == 0
        return sediment_and_flow_thickness_total_decompacted
    end
    if iuse_downhill_diffusion == 1 || use_magma_flush
        @timeit_memit "Finished updating topography using sediment transport" begin
            (
                sediment_and_flow_thickness_total_decompacted
            ) = run_sediment_transport_model!(
                model, compaction_correction_type="variable_property")
        end
    end
    return sediment_and_flow_thickness_total_decompacted
end

function update_topography_using_lava_flow!(
    model::ModelData,
    sediment_and_flow_thickness_total_decompacted::Union{Vector{Float64}, Nothing}
)::Nothing
    iuse_downhill_diffusion = model.topography.parameters.model_options.iuse_downhill_diffusion.value
    use_magma_flush = CheckMagmaFlush.check_for_initial_magma_flush(model)
    if get_iuse_topo(model) == 0
        return nothing
    end
    @timeit_memit "Finished updating topography using lava flow" begin
        if iuse_downhill_diffusion == 1 || use_magma_flush
            lava_flow_manager(model, sediment_and_flow_thickness_total_decompacted, use_magma_flush)
        else
            # For this option erosion and sedimentation was taken into account by
            # prescribed advection.
            lava_flow_manager(model, nothing, use_magma_flush)
        end
    end
end

""" Apply marker transformations for updated topography.

Apply surface processes for erosion, sedimentation, water/air boundary
and extrusion using water level and topography marker chain. Extrusion
thickness is used to define a layer of sediment with a layer of
lava flow material on top of it.
"""
function transform_markers_for_surface_processes!(model::ModelData)::Nothing
    if get_iuse_topo(model) == 1
        @timeit_memit "Finished transforming markers for surface processes" begin
            apply_marker_transformations!(model)
            zero_extrusion_thickness!(model)
        end
    end
    return nothing
end

""" Zero extrusion thickness at topography grid nodes.

Extrusion thickness is zeroed out after marker composition is updated for
erosion, sedimentation and extrusion.
"""
function zero_extrusion_thickness!(model::ModelData)::Nothing
    gridt = model.topography.arrays.gridt.array
    toponum = size(gridt, 2)
    for j in 1:toponum
        gridt[7, j] = 0.0
    end
    return nothing
end

""" Apply surface processes model by transforming markers.

Surface processes include base level, sedimentation, erosion and
volcanic extrusion.

The time the process occurred assigned to the marker as a marker age. Note
this age is actually marker time since the start of the model run. The
name age should be updated to model_time to avoid confusion.

When this function is called topography has been updated for sedimentation,
erosion and volcanic extrusion. Markers are transformed based on the
relationship between updated topography marker chain and marker material
type (matid).

This function assumes that volcanics are extruded after sedimentation for
a given time step. The calculated extrusion thickness is used to determine
the domain of the volcanic material.
"""
function apply_marker_transformations!(
    model::ModelData,
)::Nothing

    iuse_extrusion = model.melting.parameters.extrusion.iuse_extrusion.value
    iuse_extrusion = update_use_extrusion_for_magma_flush(model, iuse_extrusion)

    age_ma = calculate_age_ma(model)

    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marker_matid = model.markers.arrays.material.marker_matid.array

    y_sealevel = model.topography.parameters.sealevel.y_sealevel.value
    dx_topo = model.topography.parameters.topo_grid.dx_topo.value
    toponum = model.topography.parameters.topo_grid.toponum.value

    gridt = model.topography.arrays.gridt.array

    erosion_flag, sedimentation_flag = get_erosion_and_sedimentation_flags(model)

    matids_sticky_air = model.materials.dicts.matid_types["StickyAir"]
    matids_sticky_water = model.materials.dicts.matid_types["StickyWater"]
    matids_sediments = model.materials.dicts.matid_types["Sediment"]
    matids_volcanics = model.materials.dicts.matid_types["ExtrudedGabbroicMagma"]

    ifind_rock_above_topo = 0
    marknum = model.markers.parameters.distribution.marknum.value
    
    Threads.@threads for imarker in 1:marknum
        @inbounds begin
            x_marker = marker_x[imarker]
            y_marker = marker_y[imarker]
            matid = marker_matid[imarker]
        end
        ixn, dx_dist = TopoInterpolation.find_closest_topography_node_to_marker(
            x_marker, gridt, dx_topo, toponum)
        y_topo = TopoInterpolation.calculate_topography(dx_dist, gridt, ixn)
        # Get extrusive thickness
        if iuse_extrusion == 1
            extrusion_thickness = TopoInterpolation.calculate_extrusion_thickness(
                dx_dist, gridt, ixn)
        else
            extrusion_thickness = 0.0
        end
        if air_in_water(matid, y_marker, y_sealevel, matids_sticky_air)
            MarkerTransformation.transform_marker_to_water!(
                model, imarker, age_ma, matids_sticky_water)
        end
        if water_in_air(matid, y_marker, y_sealevel, matids_sticky_water)
            MarkerTransformation.transform_marker_to_air!(
                model, imarker, age_ma, matids_sticky_air)
        end
        if air_in_sediment(matid, y_marker, y_topo, matids_sticky_air)
            MarkerTransformation.transform_marker_to_sediment_or_volcanics!(
                model, imarker, age_ma, y_topo, y_marker, iuse_extrusion,
                extrusion_thickness, sedimentation_flag, matids_sediments,
                matids_volcanics
            )
        end
        if water_in_sediment(matid, y_marker, y_topo, matids_sticky_water)
            MarkerTransformation.transform_marker_to_sediment_or_volcanics!(
                model, imarker, age_ma, y_topo, y_marker, iuse_extrusion,
                extrusion_thickness, sedimentation_flag, matids_sediments,
                matids_volcanics
            )
        end
        if erosion(matid, y_marker, y_topo, matids_sticky_air, matids_sticky_water)
            if erosion_flag == 0
                if ifind_rock_above_topo == 0
                    println(
                        "!!! WARNING !!! Erosion is turned off but rock " *
                        "markers are located above topography marker chain. " *
                        "This may indicate instabilities in sticky layer. " *
                        "Erosion is allowed to occur to stabilize the sticky " *
                        "layer."
                    )
                    ifind_rock_above_topo = 1
                    erosion_flag = 1
                end
            end
            if erosion_flag == 1
                MarkerTransformation.transform_rocks_to_water_or_air!(
                    model, y_sealevel, imarker, age_ma,
                    matids_sticky_air, matids_sticky_water)
            end
        end
    end
end

""" Get erosion and sedimentation flags.

If flags are equal to zero markers will not be transformed due to
sedimentation or erosion.

Returns
-------
erosion_flag: Int
    Flag for erosion: 0 = no erosion, 1 = erosion.

sedimentation_flag: Int
    Flag for sedimentation: 0 = no sedimentation, 1 = sedimentation.
"""
@inline function get_erosion_and_sedimentation_flags(
    model::ModelData
)::Tuple{Int, Int}
    erosion_rate = model.topography.parameters.depo_and_erosion_rates.erosion_rate.value
    sedimentation_rate = model.topography.parameters.depo_and_erosion_rates.sedimentation_rate.value
    iuse_downhill_diffusion = model.topography.parameters.model_options.iuse_downhill_diffusion.value

    erosion_flag = 0
    if erosion_rate > 0.0 || iuse_downhill_diffusion == 1
        erosion_flag = 1
    end

    sedimentation_flag = 0
    if sedimentation_rate > 0.0 || iuse_downhill_diffusion == 1
        sedimentation_flag = 1
    end

    return erosion_flag, sedimentation_flag
end

@inline function air_in_water(
    matid::Int16,
    y_marker::Float64,
    y_sealevel::Float64,
    matids_sticky_air::Vector{Int16}
)::Bool
    check = false
    if matid in matids_sticky_air && y_marker >= y_sealevel
        check = true
    end
    return check
end

@inline function water_in_air(
    matid::Int16,
    y_marker::Float64,
    y_sealevel::Float64,
    matids_sticky_water::Vector{Int16}
)::Bool
    check = false
    if matid in matids_sticky_water && y_marker < y_sealevel
        check = true
    end
    return check
end

@inline function air_in_sediment(
    matid::Int16,
    y_marker::Float64,
    y_topo::Float64,
    matids_sticky_air::Vector{Int16}
)::Bool
    check = false
    if matid in matids_sticky_air && y_marker > y_topo
        check = true
    end
    return check
end

@inline function water_in_sediment(
    matid::Int16,
    y_marker::Float64,
    y_topo::Float64,
    matids_sticky_water::Vector{Int16}
)::Bool
    check = false
    if matid in matids_sticky_water && y_marker > y_topo
        check = true
    end
    return check
end

@inline function erosion(
    matid::Int16,
    y_marker::Float64,
    y_topo::Float64,
    matids_sticky_air::Vector{Int16},
    matids_sticky_water::Vector{Int16}
)::Bool
    sticky = false
    if matid in matids_sticky_air || matid in matids_sticky_water
        sticky = true
    end
    check = false
    if !sticky && y_marker < y_topo
        check = true
    end
    return check
end

@inline function calculate_age_ma(model::ModelData)::Float64
    timesum = model.timestep.parameters.main_time_loop.timesum.value
    age_ma = timesum/model.conversion.parameters.sec_per_Myr.value
    return age_ma
end

""" Reset marker compaction properties.

Reset marker initial porosity ``marker_porosity_initial`` and marker
porosity decay depth ``marker_decay_depth`` after markers
are transformed and recycled.
"""
function reset_marker_compaction_properties!(
    model::ModelData
)::Nothing
    @timeit_memit "Finished resetting marker compaction properties" begin
        reset_compaction_properties!(model)
    end
    return nothing
end

function apply_salt_deposition_model!(
    model::ModelData
)::Nothing
    SaltDeposition.apply_salt_deposition!(model)
    return nothing
end


end # module 