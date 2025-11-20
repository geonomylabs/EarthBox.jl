module VcycleGroup

import EarthBox.Parameters: ParameterInt, ParameterFloat
import EarthBox.Parameters: print_parameter
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

mutable struct Vcycle <: AbstractParameterGroup
    levelnum::ParameterInt
    nvcycles::ParameterInt
    convergence_criterion::ParameterFloat
    max_levels::ParameterInt
    ivcycle::ParameterInt
    iglobal_update::ParameterInt
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}

    function Vcycle(
        levelnum::ParameterInt,
        nvcycles::ParameterInt,
        convergence_criterion::ParameterFloat,
        max_levels::ParameterInt,
        ivcycle::ParameterInt,
        iglobal_update::ParameterInt
    )
        obj_list = Union{ParameterFloat, ParameterInt}[]
        new(levelnum, nvcycles, convergence_criterion, max_levels, ivcycle, iglobal_update, obj_list)
    end
end

function Vcycle(nvcycles::Int64, levelnum::Int64, max_levels::Int64)::Vcycle
    data = Vcycle(
        ParameterInt(levelnum, "levelnum", "None", "Number of multigrid resolution levels"),
        ParameterInt(nvcycles, "nvcycles", "None", "Total number of multigrid iteration cycles"),
        ParameterFloat(1e-8, "convergence_criterion", "None", "Convergence criterion."),
        ParameterInt(max_levels, "max_levels", "None", "Maximum number of multigrid levels"),
        ParameterInt(1, "ivcycle", "None", "Number of multigrid iteration cycles"),
        ParameterInt(1, "iglobal_update", "None", "Number of multigrid iteration cycles multi-multigrid iteration")
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end