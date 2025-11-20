module Options

import EarthBox.OptionTools: OptionState
import EarthBox.OptionTools: make_option_names

const option_ids = Dict{Symbol, Int}(
    :ViscoelasticStressForecast => 0,
    :PureElasticStressForecast => 1
)

const option_names = make_option_names(option_ids)

function get_options()::Dict{Int, OptionState}
    Dict{Int, OptionState}(
        option_ids[option_names.ViscoelasticStressForecast] => 
            OptionState(
                option_name=string(option_names.ViscoelasticStressForecast),
                description=(
                    "Viscoelastic plasticity model where stress is forecasted "
                    *"using a viscoelastic stress forecast."
                    ),
            ),
        option_ids[option_names.PureElasticStressForecast] => 
            OptionState(
                option_name=string(option_names.PureElasticStressForecast),
                description=(
                    "Viscoelastic plasticity model where stress is forecasted "
                    *"using a purely elastic stress forecast."
                    ),
            )
    )
end

end # module 