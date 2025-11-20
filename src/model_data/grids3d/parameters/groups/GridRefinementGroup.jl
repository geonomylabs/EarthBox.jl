module GridRefinementGroup

import EarthBox.Parameters: ParameterFloat, ParameterInt
import EarthBox.Parameters: print_parameter
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

mutable struct GridRefinement <: AbstractParameterGroup
    dx_highres::ParameterFloat
    xo_highres::ParameterFloat
    ixo_highres::ParameterInt
    xf_highres::ParameterFloat
    dz_highres::ParameterFloat
    zo_highres::ParameterFloat
    izo_highres::ParameterInt
    zf_highres::ParameterFloat
    dy_highres::ParameterFloat
    yf_highres::ParameterFloat
    iuse_trench::ParameterInt
    trench_location::ParameterFloat
    iuse_refinement_delay::ParameterInt
    refinement_time::ParameterFloat
    refinement_flag::ParameterInt
    iuse_refinement_gap::ParameterInt
    refinement_gap_start_time::ParameterFloat
    refinement_gap_end_time::ParameterFloat
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}

    function GridRefinement(
        dx_highres::ParameterFloat,
        xo_highres::ParameterFloat,
        ixo_highres::ParameterInt,
        xf_highres::ParameterFloat,
        dz_highres::ParameterFloat,
        zo_highres::ParameterFloat,
        izo_highres::ParameterInt,
        zf_highres::ParameterFloat,
        dy_highres::ParameterFloat,
        yf_highres::ParameterFloat,
        iuse_trench::ParameterInt,
        trench_location::ParameterFloat,
        iuse_refinement_delay::ParameterInt,
        refinement_time::ParameterFloat,
        refinement_flag::ParameterInt,
        iuse_refinement_gap::ParameterInt,
        refinement_gap_start_time::ParameterFloat,
        refinement_gap_end_time::ParameterFloat
    )
        obj_list = Union{ParameterFloat, ParameterInt}[]
        new(
            dx_highres, xo_highres, ixo_highres, xf_highres,
            dz_highres, zo_highres, izo_highres, zf_highres,
            dy_highres, yf_highres, iuse_trench, trench_location,
            iuse_refinement_delay, refinement_time, refinement_flag,
            iuse_refinement_gap, refinement_gap_start_time,
            refinement_gap_end_time, obj_list
        )
    end
end

function GridRefinement()::GridRefinement
    data = GridRefinement(
        ParameterFloat(
            0.0, "dx_highres", "m",
            "Horizontal grid resolution in x-direction in high-resolution area."
        ),
        ParameterFloat(
            0.0, "xo_highres", "m",
            "x-location of first node of high-resolution area."
        ),
        ParameterInt(
            0, "ixo_highres", "None",
            "x-index of first node of high-resolution area."
        ),
        ParameterFloat(
            0.0, "xf_highres", "m",
            "x-location of last node of high-resolution area"
        ),
        ParameterFloat(
            0.0, "dz_highres", "m",
            "Horizontal grid resolution in z-direction in high-resolution area."
        ),
        ParameterFloat(
            0.0, "zo_highres", "m",
            "z-location of first node of high-resolution area."
        ),
        ParameterInt(
            0, "izo_highres", "None",
            "z-index of first node of high-resolution area."
        ),
        ParameterFloat(
            0.0, "zf_highres", "m",
            "z-location of last node of high-resolution area."
        ),
        ParameterFloat(
            0.0, "dy_highres", "m",
            "Vertical grid resolution in high-resolution area."
        ),
        ParameterFloat(
            0.0, "yf_highres", "m",
            "y-location of final node of high-resolution area."
        ),
        ParameterInt(
            0, "iuse_trench", "None",
            "Flag to use trench to update high-resolution domain."
        ),
        ParameterFloat(
            0.0, "trench_location", "m",
            "Location of trench used to update high-resolution domain."
        ),
        ParameterInt(
            0, "iuse_refinement_delay", "None",
            "Flag to use refinement delay."
        ),
        ParameterFloat(
            0.0, "refinement_time", "Myr",
            "Time to start refinement."
        ),
        ParameterInt(
            0, "refinement_flag", "None",
            "Flag indicating that refinement has been applied: 0 = not applied, 1 = applied"
        ),
        ParameterInt(
            0, "iuse_refinement_gap", "None",
            "Flag to use refinement gap."
        ),
        ParameterFloat(
            0.0, "refinement_gap_start_time", "Myr",
            "Start time of refinement gap."
        ),
        ParameterFloat(
            4000.0, "refinement_gap_end_time", "Myr",
            "End time of refinement gap."
        )
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module 