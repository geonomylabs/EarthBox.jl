module ParameterCalculator

import EarthBox.PrintFuncs: print_info
import EarthBox.DictUtils: check_for_missing_key_names
import EarthBox.EarthBoxDtypes: ParametersDictType
import EarthBox: TtypeCalculator
import ..Ttype.TtypeGridData: TtypeGridDataState

function calculate_number_of_basic_nodes_x(ttype_grid_data::TtypeGridDataState)::Int64
    return TtypeCalculator.calculate_number_of_basic_nodes_x(
        ttype_grid_data.xsize,
        ttype_grid_data.xo_highres,
        ttype_grid_data.xf_highres,
        ttype_grid_data.dx_highres,
        ttype_grid_data.dx_lowres)
end

function calculate_number_of_basic_nodes_y(ttype_grid_data::TtypeGridDataState)::Int64
    return TtypeCalculator.calculate_number_of_basic_nodes_y(
        ttype_grid_data.ysize,
        ttype_grid_data.yf_highres,
        ttype_grid_data.dy_highres,
        ttype_grid_data.dy_lowres)
end

function calculate_number_of_markers_per_cell_x(ttype_grid_data::TtypeGridDataState)::Float64
    return TtypeCalculator.calculate_number_of_markers_per_cell_x(
        ttype_grid_data.xsize,
        ttype_grid_data.xo_highres,
        ttype_grid_data.xf_highres,
        ttype_grid_data.dx_highres,
        ttype_grid_data.dx_lowres,
        ttype_grid_data.dx_marker)
end

function calculate_number_of_markers_per_cell_y(ttype_grid_data::TtypeGridDataState)::Float64
    return TtypeCalculator.calculate_number_of_markers_per_cell_y(
        ttype_grid_data.ysize,
        ttype_grid_data.yf_highres,
        ttype_grid_data.dy_highres,
        ttype_grid_data.dy_lowres,
        ttype_grid_data.dy_marker)
end

end # module

