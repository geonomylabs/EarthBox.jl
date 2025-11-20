module Volcanism

import EarthBox.PrintFuncs: print_info
import EarthBox.ModelDataContainer: ModelData
import ..MeltCheck: melt_available

""" Update number of magma markers.

# Returns
- nmarkers_magma::Int
    - Number of magma markers.
"""
@inline function update_number_of_magma_markers_for_extrusion(
    nmarkers_magma::Int,
    nmarkers_volcanics::Int,
    iuse_extrusion::Int64
)::Int
    if iuse_extrusion == 1
        nmarkers_magma = nmarkers_magma - nmarkers_volcanics
    end
    return nmarkers_magma
end

""" Calculate number of volcanic markers.

# Returns
- nmarkers_volcanics::Int
    - Number of volcanic markers.
"""
function calculate_number_of_volcanic_markers(
    model::ModelData,
    nmarkers_magma::Int,
    characteristic_magmatic_crust_height::Float64,
    mantle_melting_mat_ids::Vector{Int16};
    xstart::Float64=-1e39,
    xend::Float64=1e39,
    initial_magma_flush_steps::Int=0,
    magma_flush_factor::Float64=0.8,
    min_extrusion_depth_meters::Float64=1e39
)::Int
    icount_above_minimum_extrusion_depth = 
        calculate_number_of_markers_above_minimum_extrusion_depth(
            model, mantle_melting_mat_ids, min_extrusion_depth_meters,
            xstart=xstart, xend=xend
        )

    extrusion_volume_factor = calculate_extrusion_volume_factor(
        model, characteristic_magmatic_crust_height
    )
    println("   >> extrusion_volume_factor: ", extrusion_volume_factor)

    ntimestep = model.timestep.parameters.main_time_loop.ntimestep.value
    if magma_flush_event(ntimestep, initial_magma_flush_steps)
        println("   >> Using initial magma flush at timestep: ", ntimestep)
        println("         magam flush factor: ", magma_flush_factor)
        extrusion_volume_factor = magma_flush_factor
    end

    if icount_above_minimum_extrusion_depth > 0
        nmarkers_volcanics = floor(Int, nmarkers_magma * extrusion_volume_factor)
    else
        nmarkers_volcanics = 0
    end

    print_info = false
    if print_info
        print_volcanics_info(
            icount_above_minimum_extrusion_depth,
            extrusion_volume_factor,
            nmarkers_magma, nmarkers_volcanics
        )
    end

    return nmarkers_volcanics
end

""" Calculate extrusion volume factor.
"""
function calculate_extrusion_volume_factor(
    model::ModelData,
    characteristic_magmatic_crust_height::Float64
)::Float64
    extrusion_volume_factor_min              = model.melting.parameters.extrusion.extrusion_volume_factor.value
    extrusion_volume_factor_max              = model.melting.parameters.extrusion.extrusion_volume_factor_max.value
    characteristic_magmatic_crust_height_min = model.melting.parameters.extrusion.characteristic_magmatic_crust_height_min.value
    characteristic_magmatic_crust_height_max = model.melting.parameters.extrusion.characteristic_magmatic_crust_height_max.value

    if characteristic_magmatic_crust_height <= characteristic_magmatic_crust_height_min
        extrusion_volume_factor = extrusion_volume_factor_min
    elseif characteristic_magmatic_crust_height >= characteristic_magmatic_crust_height_max
        extrusion_volume_factor = extrusion_volume_factor_max
    else
        extrusion_volume_factor = (
            extrusion_volume_factor_min
            + (extrusion_volume_factor_max - extrusion_volume_factor_min)
            / (characteristic_magmatic_crust_height_max - characteristic_magmatic_crust_height_min)
            * (characteristic_magmatic_crust_height - characteristic_magmatic_crust_height_min)
        )
    end
    return extrusion_volume_factor
end

""" Check if magma flush event should occur.
"""
function magma_flush_event(
    ntimestep::Int,
    initial_magma_flush_steps::Int
)::Bool
    check = false
    if initial_magma_flush_steps > 0
        if ntimestep <= initial_magma_flush_steps
            check = true
        end
    end
    return check
end

""" Print information about volcanic markers.
"""
function print_volcanics_info(
    icount_above_minimum_extrusion_depth::Int,
    volcanic_volume_fraction::Float64,
    nmarkers_magma::Int,
    nmarkers_volcanics::Int
)::Nothing
    print_info("", level=1)
    print_info(">> icount_above_minimum_extrusion_depth: $(icount_above_minimum_extrusion_depth)", level=1)
    print_info(">> volcanic_volume_fraction: $(volcanic_volume_fraction)", level=1)
    print_info(">> nmarkers_magma: $(nmarkers_magma)", level=1)
    print_info(">> nmarkers_volcanics: $(nmarkers_volcanics)", level=1)
    print_info("", level=1)
    return nothing
end

""" Calculate number of markers available for extraction above min depth.

If this number is zero, no volcanic markers will be created. This
effectively activates volcanism only when there are molten markers at a
shallow enough depth.

# Returns
- icount_below_minimum_extrusion_depth::Int
    - Number of markers available for melt extraction that are above the
    minimum extraction depth.
"""
function calculate_number_of_markers_above_minimum_extrusion_depth(
    model::ModelData,
    mantle_melting_mat_ids::Vector{Int16},
    min_extrusion_depth_meters::Float64;
    xstart::Float64=-1e39,
    xend::Float64=1e39
)::Int
    marker_x                    = model.markers.arrays.location.marker_x.array
    marker_y                    = model.markers.arrays.location.marker_y.array
    marker_matid                = model.markers.arrays.material.marker_matid.array
    marker_meltfrac             = model.markers.arrays.melt.marker_meltfrac.array
    marker_extractable_meltfrac = model.markers.arrays.melt.marker_extractable_meltfrac.array

    marknum = model.markers.parameters.distribution.marknum.value
    icount_above_minimum_extrusion_depth = 0
    for imarker in 1:marknum
        @inbounds begin
            x_marker = marker_x[imarker]
        end
        if xstart <= x_marker <= xend
            @inbounds matid = marker_matid[imarker]
            if matid in mantle_melting_mat_ids
                @inbounds begin
                    meltfrac = marker_meltfrac[imarker]
                    extractable_meltfrac = marker_extractable_meltfrac[imarker]
                end
                if melt_available(meltfrac, extractable_meltfrac)
                    @inbounds y_marker = marker_y[imarker]
                    if y_marker <= min_extrusion_depth_meters
                        icount_above_minimum_extrusion_depth += 1
                    end
                end
            end
        end
    end
    return icount_above_minimum_extrusion_depth
end

""" Calculate volcanic extrusion volume

# Returns
- extrusion_volume::Float64
    - Volume of extruded material in m^3.
"""
function calculate_extrusion_volume(
    model::ModelData,
    nmarkers_volcanics::Int,
    iuse_extrusion::Int64
)::Float64
    mxstep = model.markers.parameters.distribution.mxstep.value
    mystep = model.markers.parameters.distribution.mystep.value
    extrusion_volume = 0.0
    if iuse_extrusion == 1
        # volume of extruded material in m^3
        extrusion_volume = nmarkers_volcanics * mxstep * mystep
    end
    return extrusion_volume
end

end # module 