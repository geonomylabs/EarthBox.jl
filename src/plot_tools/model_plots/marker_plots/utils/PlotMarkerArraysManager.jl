module PlotMarkerArraysManager

mutable struct PlotMarkerArrays
    marker_x_km::Vector{Float64}
    marker_y_km::Vector{Float64}
    marker_matid::Vector{Float64}
    marker_strain_plastic::Vector{Float64}
    marker_strain_rate_plastic::Vector{Float64}
    marker_meltfrac::Vector{Float64}
    marker_extractable_meltfrac::Vector{Float64}
    marker_extracted_meltfrac::Vector{Float64}
    marker_pfailure::Vector{Float64}
    marker_TK::Vector{Float64}
    marker_pr::Vector{Float64}
    marker_age::Vector{Float64}
    marker_serpentinization::Vector{Float64}
    marker_rho::Vector{Float64}
end

function PlotMarkerArrays()
    return PlotMarkerArrays(
        zeros(Float64, 1),
        zeros(Float64, 1),
        zeros(Float64, 1),
        zeros(Float64, 1),
        zeros(Float64, 1),
        zeros(Float64, 1),
        zeros(Float64, 1),
        zeros(Float64, 1),
        zeros(Float64, 1),
        zeros(Float64, 1),
        zeros(Float64, 1),
        zeros(Float64, 1),
        zeros(Float64, 1),
        zeros(Float64, 1)
    )
end

function set_marker_arrays!(
    pma_state::PlotMarkerArrays,
    marker_data_dict::Dict{String, Vector{Float64}}
)::Nothing
    pma_state.marker_x_km = marker_data_dict["marker_x"]
    pma_state.marker_y_km = marker_data_dict["marker_y"]
    pma_state.marker_matid = marker_data_dict["marker_matid"]
    pma_state.marker_strain_plastic = marker_data_dict["marker_strain_plastic"]
    pma_state.marker_strain_rate_plastic = marker_data_dict["marker_strain_rate_plastic"]
    pma_state.marker_meltfrac = marker_data_dict["marker_meltfrac"]
    pma_state.marker_extractable_meltfrac = marker_data_dict["marker_extractable_meltfrac"]
    pma_state.marker_extracted_meltfrac = marker_data_dict["marker_extracted_meltfrac"]
    pma_state.marker_pfailure = marker_data_dict["marker_pfailure"]
    pma_state.marker_TK = marker_data_dict["marker_TK"]
    pma_state.marker_pr = marker_data_dict["marker_pr"]
    pma_state.marker_age = marker_data_dict["marker_age"]
    pma_state.marker_serpentinization = marker_data_dict["marker_serpentinization"]
    pma_state.marker_rho = marker_data_dict["marker_rho"]
    return nothing
end

end # module
