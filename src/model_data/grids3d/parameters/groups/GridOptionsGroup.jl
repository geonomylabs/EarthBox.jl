module GridOptionsGroup

import EarthBox.Parameters: ParameterInt, ParameterStr
import EarthBox.Parameters: print_parameter
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

mutable struct GridOptions <: AbstractParameterGroup
    itype_grid::ParameterInt
    stype_grid::ParameterStr
    obj_list::Vector{ParameterInt}

    function GridOptions(
        itype_grid::ParameterInt,
        stype_grid::ParameterStr
    )
        obj_list = ParameterInt[]
        new(itype_grid, stype_grid, obj_list)
    end
end

function GridOptions()::GridOptions
    data = GridOptions(
        ParameterInt(
            0, "itype_grid", "None",
            "Grid options: 0 = uniform; 1, 2, 3 = non-uniform"
        ),
        ParameterStr(
            "None", "stype_grid", "None",
            "Grid option name"
        )
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module 