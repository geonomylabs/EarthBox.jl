module Options

import EarthBox.OptionTools: OptionState
import EarthBox.OptionTools: make_option_names

const option_ids = Dict{Symbol, Int}(
    :Liao14 => 0,
    :Blankenbach89 => 1,
    :Naliboff17 => 2
)

const option_names = make_option_names(option_ids)

function get_options()::Dict{Int, OptionState}
    Dict{Int, OptionState}(
        option_ids[option_names.Liao14] =>
            OptionState(
                option_name=string(option_names.Liao14),
                description="This density model type is a function of "
                *"temperature and pressure following the equation described "
                *"by Liao et al. (2014)."
            ),
        option_ids[option_names.Blankenbach89] =>
            OptionState(
                option_name=string(option_names.Blankenbach89),
                description="This density model type is a function of temperature "
                *"following the equation described by Blankenbach (1989)."
            ),
        option_ids[option_names.Naliboff17] =>
            OptionState(
                option_name=string(option_names.Naliboff17),
                description="This density model uses 600 degrees Celsius as a "
                *"reference and does not include pressure."
            )
    )
end

end # module 