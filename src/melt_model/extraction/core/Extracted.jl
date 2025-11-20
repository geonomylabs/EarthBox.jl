module Extracted

import EarthBox.ModelDataContainer: ModelData
import ..MeltCheck: melt_available

""" Update total extracted melt fraction for markers.

The extracted melt fraction Mextracted_i integrates incremental amounts
of melt that are extracted at each time step.

First the incremental exactable melt fraction Mextractable_i is calculated
during a previous step using the following equation:

    Mextractable_i = Mfraction_i - Mextracted_i-1                   (eq. 1)

where Mfraction_i is the total melt fraction at the current time step i and
Mextracted_i-1 is the extracted melt fraction from the previous time step
i. The the updated extracted melt fraction, which is the focus of this
function, is calculated using the following equation:

    Mextracted_i = Mextracted_i-1 + Mextractable_i                  (eq. 2)

where Mextracted_i-1 is the previous extracted melt fraction and
Mextractable_i is the current extractable melt fraction calculated with
eq. 1.

# Updated Array
- model.markers.arrays.melt.marker_extracted_meltfrac.array::Vector{Float64}
    - Extracted melt fraction of marker.
"""
function update_extracted_meltfrac!(
    model::ModelData,
    mantle_melting_mat_ids::Vector{Int16}
)::Nothing
    marker_matid = model.markers.arrays.material.marker_matid.array

    melt_arrays = model.markers.arrays.melt
    marker_meltfrac = melt_arrays.marker_meltfrac.array
    marker_extracted_meltfrac = melt_arrays.marker_extracted_meltfrac.array
    marker_extractable_meltfrac = melt_arrays.marker_extractable_meltfrac.array

    marknum = model.markers.parameters.distribution.marknum.value
    
    for imarker in 1:marknum
        if marker_matid[imarker] in mantle_melting_mat_ids
            meltfrac = marker_meltfrac[imarker]
            extractable_meltfrac = marker_extractable_meltfrac[imarker]
            if melt_available(meltfrac, extractable_meltfrac)
                marker_extracted_meltfrac[imarker] += extractable_meltfrac
            end
        end
    end
    return nothing
end

end # module 