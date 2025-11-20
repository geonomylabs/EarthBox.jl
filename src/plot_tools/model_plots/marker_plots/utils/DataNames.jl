module DataNames

import EarthBox.JLDTools: get_jld_marker_data

Base.@kwdef struct MarkerDataNames
    marker_x::String = "marker_x"
    marker_y::String = "marker_y"
    marker_matid::String = "marker_matid"
    marker_plastic_strain::String = "marker_strain_plastic"
    marker_strain_rate_plastic::String = "marker_strain_rate_plastic"
    marker_meltfrac::String = "marker_meltfrac"
    marker_extractable_meltfrac::String = "marker_extractable_meltfrac"
    marker_extracted_meltfrac::String = "marker_extracted_meltfrac"
    marker_pfailure::String = "marker_pfailure"
    marker_TK::String = "marker_TK"
    marker_pr::String = "marker_pr"
    marker_age::String = "marker_age"
    marker_serpentinization::String = "marker_serpentinization"
    marker_rho::String = "marker_rho"
end

function get_list(mdn::MarkerDataNames)
    return [
        mdn.marker_x,
        mdn.marker_y,
        mdn.marker_matid,
        mdn.marker_plastic_strain,
        mdn.marker_strain_rate_plastic,
        mdn.marker_meltfrac,
        mdn.marker_extractable_meltfrac,
        mdn.marker_extracted_meltfrac,
        mdn.marker_pfailure,
        mdn.marker_TK,
        mdn.marker_pr,
        mdn.marker_age,
        mdn.marker_serpentinization,
        mdn.marker_rho
    ]
end

"""
    get_data(mdn::MarkerDataNames, itime_step::Int, dir_path::String, 
             dataname::String)

Get marker data from JLD file for the specified time step and data name.
Returns a tuple of (time, data_array).
"""
function get_data(
    mdn::MarkerDataNames, 
    itime_step::Int, 
    dir_path::String, 
    dataname::String
)
    check_dataname(mdn, dataname)
    return get_jld_marker_data(itime_step, dir_path, dataname)
end

"""
    check_dataname(mdn::MarkerDataNames, dataname::String)

Check if the provided data name is valid. Raises an error if invalid.
"""
function check_dataname(mdn::MarkerDataNames, dataname::String)
    valid_datanames = get_list(mdn)
    if !(dataname in valid_datanames)
        println(">> Valid data names:")
        for name in valid_datanames
            println("    $name")
        end
        error("Invalid data name: $dataname. See list of valid names above.")
    end
end

end # module