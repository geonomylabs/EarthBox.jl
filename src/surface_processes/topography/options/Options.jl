module Options

import EarthBox.OptionTools: OptionState
import EarthBox.OptionTools: make_option_names
import EarthBox.ModelDataContainer: ModelData

const option_ids = Dict{Symbol, Int}(
    :RungeKuttaWithInterp => 0,
    :Antidiffusion => 1
)

const option_names = make_option_names(option_ids)

function get_options()::Dict{Int, OptionState}
    Dict{Int, OptionState}(
        option_ids[option_names.RungeKuttaWithInterp] => 
            OptionState(
                option_name=string(option_names.RungeKuttaWithInterp),
                description="Topography nodes are advected using a 4th order " *
                    "Runge Kutta scheme and the elevation of the advected " *
                    "marker chain is interpolated to Eulerian topography " *
                    "grid points.",
                bools=Dict{Symbol, Bool}()
            ),
        option_ids[option_names.Antidiffusion] => 
            OptionState(
                option_name=string(option_names.Antidiffusion),
                description="Elevation is calculated at topography nodes " *
                    "using an antidiffusion scheme.",
                bools=Dict{Symbol, Bool}()
            )
    )
end

end # module 