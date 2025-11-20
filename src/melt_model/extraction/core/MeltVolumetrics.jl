module MeltVolumetrics

import EarthBox.ModelDataContainer: ModelData
import EarthBox.ConversionFuncs: seconds_to_years
import ..MeltCheck: melt_available

""" Calculate incremental total melt volume.

The incremental melt volume is calculated by summing incremental
extractable melt fractions for each marker. Extractable melt fractions are
the difference between the current total melt fraction at a given pressure
and temperature less the extracted melt fraction. The incremental melt
volume includes residual incremental melt volume from the previous time
step.

# Arguments
- model::ModelData
    - Model data.
- mantle_melting_mat_ids::Vector{Int16}
    - Material IDs for mantle melting.
- xstart::Float64
    - Start of the domain in the x-direction (meters).
- xend::Float64
    - End of the domain in the x-direction (meters).

# Returns
- melt_volume::Float64
    - Total melt volume in terms of marker units. A volume of 1 equals the
    volume of the marker in marker units. To convert to m^3 multiply
    by the volume occupied by each marker based on average marker spacing.
"""
function calculate_melt_volume(
    model::ModelData,
    mantle_melting_mat_ids::Vector{Int16};
    xstart::Float64=-1e39,
    xend::Float64=1e39
)::Float64
    marknum = model.markers.parameters.distribution.marknum.value

    marker_matid = model.markers.arrays.material.marker_matid.array
    marker_meltfrac = model.markers.arrays.melt.marker_meltfrac.array
    marker_extractable_meltfrac = 
        model.markers.arrays.melt.marker_extractable_meltfrac.array
    marker_x = model.markers.arrays.location.marker_x.array

    melt_volume = 0.0
    for imarker in 1:marknum
        x_marker = marker_x[imarker]
        if xstart <= x_marker <= xend
            if marker_matid[imarker] in mantle_melting_mat_ids
                meltfrac = marker_meltfrac[imarker]
                extractable_meltfrac = marker_extractable_meltfrac[imarker]
                if melt_available(meltfrac, extractable_meltfrac)
                    melt_volume += extractable_meltfrac
                end
            end
        end
    end
    return melt_volume
end

""" Calculate residual melt volume.

The residual melt volume is the fractional part that cannot be translated
into a whole number of markers. This residual is added to the melt value in
the next time step.

# Arguments
- melt_volume::Float64
    - Total melt volume in marker volume units where 1 is the volume of a
    single marker.
- nmarkers_magma::Int
    - Melt volume converted to an integer thus eliminating the fractional
    component.

# Returns
- melt_volume_residual_initial::Float64
    - Residual melt volume in marker units where 1 = volume of a marker.
"""
function calculate_residual_melt_volume(
    melt_volume::Float64,
    nmarkers_magma::Int
)::Float64
    melt_volume_residual = melt_volume - nmarkers_magma

    if melt_volume_residual > 1.0
        melt_volume_residual = melt_volume_residual - floor(Int, melt_volume_residual)
    end

    return melt_volume_residual
end

""" Calculate residual melt volume.

Remove non-fractional part of melt_volume_residual_initial since this was
added to number of partial melt markers.

# Returns
- melt_volume_residual::Float64
    - Residual melt volume in marker units where 1 = volume of a marker.
"""
function correct_residual_melt_volume(
    melt_volume_residual_initial::Float64
)::Float64
    if melt_volume_residual_initial > 1.0
        melt_volume_residual = (
            melt_volume_residual_initial
            - floor(Int, melt_volume_residual_initial)
        )
    else
        melt_volume_residual = melt_volume_residual_initial
    end

    return melt_volume_residual
end

""" Update number of partially molten markers.

Non-fractional part of melt_volume_residual_initial is added to
nmarkers_partial_melt.

# Arguments
- melt_volume_residual_initial::Float64
    - Residual melt volume in marker units where 1 = volume of a marker.
- nmarkers_partial_melt::Int
    - Number of markers with partial melt.

# Returns
- nmarkers_partial_melt::Int
    - Number of markers with partial melt.
"""
function correct_number_of_partially_molten_markers(
    melt_volume_residual_initial::Float64,
    nmarkers_partial_melt::Int
)::Int
    if melt_volume_residual_initial > 1
        nmarkers_partial_melt += floor(Int, melt_volume_residual_initial)
    end
    return nmarkers_partial_melt
end

""" Add residual melt volume to total volume.

# Returns
- melt_volume::Float64
    - Total melt volume in marker volume units where 1 is the volume of a
    single marker.
"""
function add_residual_volume(
    model::ModelData,
    melt_volume::Float64
)::Float64
    melt_volume_residual = model.melting.parameters.extraction.melt_residual.value
    melt_volume += melt_volume_residual
    return melt_volume
end

""" Calculate volcanic extrusion volume

# Arguments
- model::ModelData
    - Model data structure.
- nmarkers_magma::Int
    - Number of magma markers.

# Returns
- magma_production_rate_m3_yr::Float64
    - Volume of produced magma material in m^3 per year.
- avg_marker_volume_m3::Float64
    - Average volume of a marker in m^3.
"""
function calculate_magma_production_rate(
    model::ModelData,
    nmarkers_magma::Int
)::Tuple{Float64, Float64}
    avg_marker_volume_m3 = calculate_avg_marker_volume(model)
    timestep_seconds = model.timestep.parameters.main_time_loop.timestep.value
    timestep_years = seconds_to_years(timestep_seconds)
    magma_volume = nmarkers_magma * avg_marker_volume_m3
    magma_production_rate_m3_yr = magma_volume / timestep_years
    return magma_production_rate_m3_yr, avg_marker_volume_m3
end

""" Calculate average marker volume.

# Arguments
- model::ModelData
    - Model data structure.

# Returns
- avg_marker_volume_m3::Float64
    - Average volume of a marker in m^3.
"""
function calculate_avg_marker_volume(
    model::ModelData
)::Float64
    avg_marker_x_spacing_m = model.markers.parameters.distribution.mxstep.value
    avg_marker_y_spacing_m = model.markers.parameters.distribution.mystep.value
    # assume that the average marker spacing in the z-direction is 1.0 m
    avg_marker_volume_m3 = avg_marker_y_spacing_m * avg_marker_x_spacing_m
    return avg_marker_volume_m3
end

end # module 