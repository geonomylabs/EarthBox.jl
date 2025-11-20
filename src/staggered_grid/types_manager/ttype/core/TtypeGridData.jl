module TtypeGridData

import EarthBox.CaseInputTools.CaseTypes: CaseType
import EarthBox.DictUtils: check_dictionary_names

"""
    TtypeGridDataState

Main data structure for managing t-type refined grid calculations.
T-type refinement involves a high-resolution area in the upper center of the
model domain with a low x-resolution areas on either side and a low 
y-resolution area below the upper high-resolution area.
"""
mutable struct TtypeGridDataState
    # Width of model in meters
    xsize::Float64
    # Height of model in meters
    ysize::Float64
    # Initial x-grid location of high-resolution area in meters
    xo_highres::Float64
    # Final x-grid location of high-resolution area in meters
    xf_highres::Float64
    # X-direction spacing in high-resolution area in meters
    dx_highres::Float64
    # X-direction spacing in low-resolution area in meters
    dx_lowres::Float64
    # Final y-grid location of high-resolution area in meters
    yf_highres::Float64
    # Y-direction spacing in high-resolution area in meters
    dy_highres::Float64
    # Y-direction spacing in low-resolution area in meters
    dy_lowres::Float64
    # X-direction spacing of markers in meters
    dx_marker::Float64
    # Y-direction spacing of markers in meters
    dy_marker::Float64
    # X-direction spacing of topography grid in meters
    dx_topo::Float64
end

function TtypeGridDataState(case_inputs_active::CaseType)
    return TtypeGridDataState(
        case_inputs_active["xsize"].value,
        case_inputs_active["ysize"].value,
        case_inputs_active["xo_highres"].value,
        case_inputs_active["xf_highres"].value,
        case_inputs_active["dx_highres"].value,
        case_inputs_active["dx_lowres"].value,
        case_inputs_active["yf_highres"].value,
        case_inputs_active["dy_highres"].value,
        case_inputs_active["dy_lowres"].value,
        case_inputs_active["dx_marker"].value,
        case_inputs_active["dy_marker"].value,
        case_inputs_active["dx_topo"].value
    )
end
function TtypeGridDataState(ttype_grid_parameters::Dict{String, Float64})
    check_dictionary_names(REQUIRED_KEYS, ttype_grid_parameters)
    return TtypeGridDataState(
        ttype_grid_parameters["xsize"],
        ttype_grid_parameters["ysize"],
        ttype_grid_parameters["xo_highres"],
        ttype_grid_parameters["xf_highres"],
        ttype_grid_parameters["dx_highres"],
        ttype_grid_parameters["dx_lowres"],
        ttype_grid_parameters["yf_highres"],
        ttype_grid_parameters["dy_highres"],
        ttype_grid_parameters["dy_lowres"],
        ttype_grid_parameters["dx_marker"],
        ttype_grid_parameters["dy_marker"],
        ttype_grid_parameters["dx_topo"]
    )
end

"""
    check_for_problems(grid::TtypeGridState)

Check for problems with T-type grid parameters.
"""
function check_for_problems(grid::TtypeGridDataState)
    if grid.xo_highres >= grid.xf_highres
        throw(ArgumentError("T-type refinement: xo_highres must be less than xf_highres."))
    end
    if grid.xo_highres > grid.xsize
        throw(ArgumentError("T-type refinement: xo_highres must be less than model width."))
    end
    if grid.xf_highres > grid.xsize
        throw(ArgumentError("T-type refinement: xf_highres must be less than model width."))
    end
    if grid.yf_highres > grid.ysize
        throw(ArgumentError("T-type refinement: yf_highres must be less than model height."))
    end
    if (grid.xf_highres - grid.xo_highres) > grid.xsize
        throw(ArgumentError("T-type refinement: xf_highres - xo_highres must be less than model width."))
    end
    if grid.dx_highres >= (grid.xf_highres - grid.xo_highres)
        throw(ArgumentError("T-type refinement: dx_highres must be less than xf_highres - xo_highres."))
    end
    if grid.dx_lowres >= (grid.xsize - (grid.xf_highres - grid.xo_highres))/2.0
        throw(ArgumentError("T-type refinement: dx_lowres must be less than (model width - (xf_highres - xo_highres))/2."))
    end
    if grid.dy_highres >= grid.ysize
        throw(ArgumentError("T-type refinement: dy_highres must be less than model height."))
    end
    if grid.dy_lowres >= (grid.ysize - grid.yf_highres)
        throw(ArgumentError("T-type refinement: dy_lowres must be less than model height - yf_highres."))
    end
    if grid.dx_highres <= 0.0
        throw(ArgumentError("T-type refinement: dx_highres must be greater than zero."))
    end
    if grid.dx_lowres <= 0.0
        throw(ArgumentError("T-type refinement: dx_lowres must be greater than zero."))
    end
    if grid.yf_highres <= 0.0
        throw(ArgumentError("T-type refinement: yf_highres must be greater than zero."))
    end
    if grid.dy_highres <= 0.0
        throw(ArgumentError("T-type refinement: dy_highres must be greater than zero."))
    end
    if grid.dy_lowres <= 0.0
        throw(ArgumentError("T-type refinement: dy_lowres must be greater than zero."))
    end
    if grid.dx_marker <= 0.0
        throw(ArgumentError("T-type refinement: dx_marker must be greater than zero."))
    end
    if grid.dy_marker <= 0.0
        throw(ArgumentError("T-type refinement: dy_marker must be greater than zero."))
    end
    if grid.dx_topo <= 0.0
        throw(ArgumentError("T-type refinement: dx_topo must be greater than zero."))
    end
end

end
