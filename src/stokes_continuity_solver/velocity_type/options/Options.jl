module Options

import EarthBox.OptionTools: OptionState
import EarthBox.OptionTools: make_option_names

const option_ids = Dict{Symbol, Int}(
    :ZeroVelocity => -1,
    :VelocityFromStokesSolver => 0,
    :KinematicSolidBodyRotation => 1
)

const option_names = make_option_names(option_ids)

function get_options()::Dict{Int, OptionState}
    Dict{Int, OptionState}(
        option_ids[option_names.ZeroVelocity] => 
            OptionState(
                option_name=string(option_names.ZeroVelocity),
                description="Velocity is set to zero for testing purposes."
            ),
        option_ids[option_names.VelocityFromStokesSolver] => 
            OptionState(
                option_name=string(option_names.VelocityFromStokesSolver),
                description="Velocity is calculated via solving the Stokes-continuity equations."
            ),
        option_ids[option_names.KinematicSolidBodyRotation] => 
            OptionState(
                option_name=string(option_names.KinematicSolidBodyRotation),
                description="Velocity is prescribed using a solid body rotation model."
            )
    )
end

end # module 