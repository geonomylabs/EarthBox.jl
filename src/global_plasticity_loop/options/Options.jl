module Options

import EarthBox.OptionTools: OptionState
import EarthBox.OptionTools: make_option_names

const option_ids = Dict{Symbol, Int}(
    :MarkerPlasticityLoop => 0,
    :NodalPlasticityLoop => 1
)

const option_names = make_option_names(option_ids)

function get_options()::Dict{Int, OptionState}
    Dict{Int, OptionState}(
        option_ids[option_names.MarkerPlasticityLoop] => 
            OptionState(
                option_name=string(option_names.MarkerPlasticityLoop),
                description="Plasticity is defined on markers.",
                bools=Dict{Symbol, Bool}(:marker_plasticity => true)
            ),
        option_ids[option_names.NodalPlasticityLoop] => 
            OptionState(
                option_name=string(option_names.NodalPlasticityLoop),
                description="Plasticity is defined on nodes.",
                bools=Dict{Symbol, Bool}(:nodal_plasticity => true)
            )
    )
end

end # module 