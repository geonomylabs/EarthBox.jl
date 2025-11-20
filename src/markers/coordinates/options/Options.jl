module Options

import EarthBox.OptionTools: OptionState
import EarthBox.OptionTools: make_option_names

const option_ids = Dict{Symbol, Int}(
    :Regular => 0,
    :Randomized => 1
)

const option_names = make_option_names(option_ids)

function get_options()::Dict{Int, OptionState}
    Dict{Int, OptionState}(
        option_ids[option_names.Regular] =>
            OptionState(
                option_name=string(option_names.Regular),
                description="Initial marker distribution is regular.",
            ),
        option_ids[option_names.Randomized] =>
            OptionState(
                option_name=string(option_names.Randomized),
                description="Initial marker distribution is randomized.",
            )
    )
end

end # module 