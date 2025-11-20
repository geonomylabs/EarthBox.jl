module Options

import EarthBox.OptionTools: OptionState
import EarthBox.OptionTools: make_option_names

const option_ids = Dict{Symbol, Int}(
    :Constant => 0,
    :LeftEdge => 1,
    :AveragePressure => 2
)

const option_names = make_option_names(option_ids)

function get_options()::Dict{Int, OptionState}
    Dict{Int, OptionState}(
        option_ids[option_names.Constant] => 
            OptionState(
                option_name=string(option_names.Constant),
                description="Sealevel is at a constant depth relative to the surface "
                *"of the model."
            ),
        option_ids[option_names.LeftEdge] => 
            OptionState(
                option_name=string(option_names.LeftEdge),
                description="Sea level is adjusted by calculating a isostatic base "
                *"level relative to the left-hand edge of sticky-rock interface."
            ),
        option_ids[option_names.AveragePressure] => 
            OptionState(
                option_name=string(option_names.AveragePressure),
                description="Sea level is adjusted by calculating a isostatic base "
                *"level relative to the average pressure along the base of "
                *"the model."
            )
    )
end

end # module 