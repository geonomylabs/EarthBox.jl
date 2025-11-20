module Options

import DataStructures: OrderedDict
import EarthBox.ModelDataContainer: ModelData
import EarthBox.OptionTools: OptionState
import EarthBox.OptionTools: make_option_names

const option_ids = Dict{Symbol, Int}(
    :Regular => 0,
    :Randomized => 1
)

const option_names = make_option_names(option_ids)

function get_options()::OrderedDict{Int, OptionState}
    return OrderedDict{Int, OptionState}(
        option_ids[option_names.Regular] =>
            OptionState(
                option_name=string(option_names.Regular),
                description="Initial marker friction is regular."
            ),
        option_ids[option_names.Randomized] =>
            OptionState(
                option_name=string(option_names.Randomized),
                description="Initial marker friction is randomized using `delta_fric_coef`."
            )
    )
end

end # module 