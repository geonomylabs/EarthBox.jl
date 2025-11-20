module Options

import EarthBox.OptionTools: OptionState
import EarthBox.OptionTools: make_option_names

const option_ids = Dict{Symbol, Int}(
    :UniformGrid => 0,
    :TtypeRefinedGrid => 1,
    :Grid51x51UniformRefinementAlongSides => 2,
    :Grid51x51SmoothedRefinementAlongSides => 3
)

const option_names = make_option_names(option_ids)

function get_options()::Dict{Int, OptionState}
    Dict{Int, OptionState}(
        option_ids[option_names.UniformGrid] => 
            OptionState(
                option_name=string(option_names.UniformGrid),
                description="Uniform basic grid.",
                bools=Dict{Symbol, Bool}(:uniform_grid => true)
            ),
        option_ids[option_names.TtypeRefinedGrid] => 
            OptionState(
                option_name=string(option_names.TtypeRefinedGrid),
                description="Basic grid has lateral refinement in the central part of "
                *"the model and vertical refinement in the upper part of the model "
                *"(i.e. T-type).",
            ),
        option_ids[option_names.Grid51x51UniformRefinementAlongSides] => 
            OptionState(
                option_name=string(option_names.Grid51x51UniformRefinementAlongSides),
                description="51x51 basic grid with high-resolution along sides. This "
                *"type is used for convection in a box benchmarks.",
            ),
        option_ids[option_names.Grid51x51SmoothedRefinementAlongSides] => 
            OptionState(
                option_name=string(option_names.Grid51x51SmoothedRefinementAlongSides),
                description="51x51 basic grid with smoothed high-resolution zones along "
                *"sides. This grid is used for convection in a box benchmarks.",
            )
    )
end

end # module 
