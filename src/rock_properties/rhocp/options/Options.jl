module Options

import EarthBox.OptionTools: OptionState
import EarthBox.OptionTools: make_option_names

const option_ids = Dict{Symbol, Int}(
    :Constant => 0,
    :TemperatureDependentWaples => 1
)

const option_names = make_option_names(option_ids)

function get_options()::Dict{Int, OptionState}
    Dict{Int, OptionState}(
        option_ids[option_names.Constant] =>
            OptionState(
                option_name=string(option_names.Constant),
                description="This rho*cp model type uses a constant heat capacity."
            ),
        option_ids[option_names.TemperatureDependentWaples] =>
            OptionState(
                option_name=string(option_names.TemperatureDependentWaples),
                description="This rho*cp model type uses a temperature-dependent heat " *
                "capacity following the approach of Waples and Waples (2004). " *
                "See Hantschel and Kauerauf (2010, pg. 117) for details."
            )
    )
end

end # module 