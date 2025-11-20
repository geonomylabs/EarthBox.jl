module Options

import EarthBox.OptionTools: OptionState
import EarthBox.OptionTools: make_option_names

const option_ids = Dict{Symbol, Int}(
    :NoMotion => 0,
    :SimpleAdvection => 1,
    :RungeKutta4thOrder => 4
)

const option_names = make_option_names(option_ids)

function get_options()::Dict{Int, OptionState}
    Dict{Int, OptionState}(
        option_ids[option_names.NoMotion] => 
            OptionState(
                option_name=string(option_names.NoMotion),
                description="Markers are fixed in space and time"
            ),
        option_ids[option_names.SimpleAdvection] => 
            OptionState(
                option_name=string(option_names.SimpleAdvection),
                description="Markers are advected in the velocity field using a 1st " *
                "order Runge Kutta scheme (simple advection)."
            ),
        option_ids[option_names.RungeKutta4thOrder] => 
            OptionState(
                option_name=string(option_names.RungeKutta4thOrder),
                description="Markers are advected in the velocity field using a 4th " *
                "order Runge Kutta scheme."
            )
    )
end

end # module 