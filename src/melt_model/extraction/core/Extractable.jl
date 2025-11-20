module Extractable

import EarthBox.ModelDataContainer: ModelData

""" Update incremental extractable melt fraction for markers.

# Updated Array
- model.markers.arrays.melt.marker_extractable_meltfrac.array::Vector{Float64}
    - Extractable melt fraction of marker.
"""
function update_extractable_meltfrac!(
    model::ModelData,
    mantle_melting_mat_ids::Vector{Int16}
)
    marker_matid = model.markers.arrays.material.marker_matid.array
    marker_meltfrac = model.markers.arrays.melt.marker_meltfrac.array
    marker_extracted_meltfrac = model.markers.arrays.melt.marker_extracted_meltfrac.array
    marker_extractable_meltfrac = model.markers.arrays.melt.marker_extractable_meltfrac.array

    marknum = model.markers.parameters.distribution.marknum.value

    for imarker in 1:marknum
        if marker_matid[imarker] in mantle_melting_mat_ids
            (
                marker_extractable_meltfrac[imarker]
            ) = calculate_incremental_extractable_melt_fraction(
                marker_meltfrac[imarker],
                marker_extracted_meltfrac[imarker]
            )
        end
    end
end

""" Calculate incremental extractable melt fraction.

Extractable melt fraction is defined by the following equation:

    Mextractable_i = Mfraction_i - Mextracted_i-1

where Mfraction_i is the current calculated total melt fraction based on
temperature and pressure for marker and Mextracted_i-1 is the total
melt fraction that has been extracted from marker calculated during the
previous time step i-1.

A negative extractable melt fraction means that melt has been removed and
temperature and pressure conditions are no longer suitable for additional
melt extraction.

# Returns
- extractable_meltfrac::Float64
    - Extractable melt fraction.
"""
function calculate_incremental_extractable_melt_fraction(
    meltfrac::Float64,
    extracted_meltfrac::Float64
)::Float64
    extractable_meltfrac = meltfrac - extracted_meltfrac
    return extractable_meltfrac
end

end # module 