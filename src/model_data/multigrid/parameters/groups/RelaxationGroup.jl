module RelaxationGroup

import EarthBox.Parameters: ParameterFloat
import EarthBox.Parameters: print_parameter
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

mutable struct Relaxation <: AbstractParameterGroup
    krelaxs::ParameterFloat
    vkoef::ParameterFloat
    krelaxc::ParameterFloat
    pkoef::ParameterFloat
    obj_list::Vector{Union{ParameterFloat}}

    function Relaxation(
        krelaxs::ParameterFloat,
        vkoef::ParameterFloat,
        krelaxc::ParameterFloat,
        pkoef::ParameterFloat
    )
        obj_list = Union{ParameterFloat}[]
        new(krelaxs, vkoef, krelaxc, pkoef, obj_list)
    end
end

function Relaxation()::Relaxation
    data = Relaxation(
        ParameterFloat(0.9, "krelaxs", "None", "Relaxation coefficient for Gauss-Seidel iteration (Stokes equations)"),
        ParameterFloat(1.0, "vkoef", "None", "Relaxation coefficient for Gauss-Seidel iteration (Stokes equations)"),
        ParameterFloat(0.3, "krelaxc", "None", "Relaxation coefficient for Gauss-Seidel iteration (Continuity equation)"),
        ParameterFloat(1.0, "pkoef", "None", "Relaxation coefficient for Gauss-Seidel iteration (Continuity equation)")
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end