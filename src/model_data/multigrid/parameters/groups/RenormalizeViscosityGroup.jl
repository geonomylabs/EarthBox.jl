module RenormalizeViscosityGroup

import EarthBox.Parameters: ParameterInt, ParameterFloat
import EarthBox.Parameters: print_parameter
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

mutable struct RenormalizeViscosity <: AbstractParameterGroup
    viscosity_min_cur::ParameterFloat
    viscosity_max_cur::ParameterFloat
    viscosity_max_start::ParameterFloat
    viscosity_min_factor::ParameterFloat
    viscosity_max_factor::ParameterFloat
    viscosity_max_factor_mm::ParameterFloat
    viscosity_min_o::ParameterFloat
    viscosity_max_o::ParameterFloat
    nvcycles_viscosity_jump::ParameterInt
    nvcycles_viscosity_jump_mm::ParameterInt
    iviscosity_jump::ParameterInt
    iviscosity_max_overstep::ParameterInt
    resmax::ParameterFloat
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}

    function RenormalizeViscosity(
        viscosity_min_cur::ParameterFloat,
        viscosity_max_cur::ParameterFloat,
        viscosity_max_start::ParameterFloat,
        viscosity_min_factor::ParameterFloat,
        viscosity_max_factor::ParameterFloat,
        viscosity_max_factor_mm::ParameterFloat,
        viscosity_min_o::ParameterFloat,
        viscosity_max_o::ParameterFloat,
        nvcycles_viscosity_jump::ParameterInt,
        nvcycles_viscosity_jump_mm::ParameterInt,
        iviscosity_jump::ParameterInt,
        iviscosity_max_overstep::ParameterInt,
        resmax::ParameterFloat
    )
        obj_list = Union{ParameterFloat, ParameterInt}[]
        new(
            viscosity_min_cur, viscosity_max_cur, viscosity_max_start, viscosity_min_factor, 
            viscosity_max_factor, viscosity_max_factor_mm, viscosity_min_o, viscosity_max_o, nvcycles_viscosity_jump, nvcycles_viscosity_jump_mm, 
            iviscosity_jump, iviscosity_max_overstep, resmax, obj_list
            )
    end
end

function RenormalizeViscosity()::RenormalizeViscosity
    data = RenormalizeViscosity(
        ParameterFloat(1e20, "viscosity_min_cur", "Pa.s", "Minimum viscosity for current iteration"),
        ParameterFloat(1e20, "viscosity_max_cur", "Pa.s", "Maximum viscosity for current iteration"),
        ParameterFloat(1e20, "viscosity_max_start", "Pa.s", "Maximum viscosity for beginning of multi-multigrid iteration"),
        ParameterFloat(1.0, "viscosity_min_factor", "None", "Factor used to decrease minimum viscosity after a finite number of steps"),
        ParameterFloat(3.3333, "viscosity_max_factor", "None", "Factor used to increase maximum viscosity after a finite number of steps"),
        ParameterFloat(10.0001, "viscosity_max_factor_mm", "None", "Factor used to increase maximum viscosity after a finite number of steps for multi-multigrid"),
        ParameterFloat(1e20, "viscosity_min_o", "Pa.s", "Minimum viscosity"),
        ParameterFloat(1e23, "viscosity_max_o", "Pa.s", "Maximum viscosity"),
        ParameterInt(15, "nvcycles_viscosity_jump", "None", "Number of cycles for viscosity change."),
        ParameterInt(3, "nvcycles_viscosity_jump_mm", "None", "Number of cycles for multi-multigrid viscosity change"),
        ParameterInt(1, "iviscosity_jump", "None", "Number of cycles for viscosity change"),
        ParameterInt(0, "iviscosity_max_overstep", "None", "Reset flag for multi-multigrid"),
        ParameterFloat(1e-6, "resmax", "None", "Maximum residual")
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end