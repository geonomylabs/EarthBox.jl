module ArrayLookup

import EarthBox.ModelDataContainer: ModelData
import EarthBox.Arrays.ArrayTypes.MarkerArrayFloat1D: MarkerArrayFloat1DState
import EarthBox.Arrays.ArrayTypes.MarkerArrayInt1D: MarkerArrayInt1DState
import ..DataNames: MarkerDataNames

mutable struct ArrayLookupStruct
    datanames::MarkerDataNames
    lookup::Dict{String, Union{MarkerArrayFloat1DState, MarkerArrayInt1DState{<:Integer}}}
end

function create_array_lookup(model::ModelData)
    datanames = MarkerDataNames()
    lookup = build_lookup(model, datanames)
    return ArrayLookupStruct(datanames, lookup)
end

function build_lookup(
    model::ModelData,
    datanames::MarkerDataNames
)::Dict{String, Union{MarkerArrayFloat1DState, MarkerArrayInt1DState{<:Integer}}}
    lookup = Dict{String, Union{MarkerArrayFloat1DState, MarkerArrayInt1DState{<:Integer}}}(
        datanames.marker_x => model.markers.arrays.location.marker_x,
        datanames.marker_y => model.markers.arrays.location.marker_y,
        datanames.marker_matid => model.markers.arrays.material.marker_matid,
        datanames.marker_plastic_strain => model.markers.arrays.strain.marker_strain_plastic,
        datanames.marker_strain_rate_plastic => model.markers.arrays.strain.marker_strain_rate_plastic,
        datanames.marker_meltfrac => model.markers.arrays.melt.marker_meltfrac,
        datanames.marker_extractable_meltfrac => model.markers.arrays.melt.marker_extractable_meltfrac,
        datanames.marker_extracted_meltfrac => model.markers.arrays.melt.marker_extracted_meltfrac,
        datanames.marker_pfailure => model.markers.arrays.rheology.marker_pfailure,
        datanames.marker_TK => model.markers.arrays.thermal.marker_TK,
        datanames.marker_pr => model.markers.arrays.pressure.marker_pr,
        datanames.marker_age => model.markers.arrays.strat.marker_age,
        datanames.marker_serpentinization => model.markers.arrays.material.marker_serpentinization,
        datanames.marker_rho => model.markers.arrays.material.marker_rho
    )
    return lookup
end

function get_array_info(
    lookup_struct::ArrayLookupStruct,
    dataname::String
)::Tuple{Vector{Float64}, String}
    marker_array = lookup_struct.lookup[dataname]
    array = marker_array.array
    units = marker_array.units
    return array, units
end

end # module
