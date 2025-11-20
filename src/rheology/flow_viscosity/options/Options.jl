module Options

import EarthBox.OptionTools: OptionState
import EarthBox.OptionTools: make_option_names

const option_ids = Dict{Symbol, Int}(
    :XCoordinateDependent => -1,
    :Isoviscous => 0,
    :Composite => 1,
    :Composite2 => 2, # This is used to maintain backward compatibility with old models
    :TemperatureDependentCouetteBenchmark => 4,
    :TemperatureDependentConvectionBenchmark => 5
)

const option_names = make_option_names(option_ids)

function get_options()::Dict{Int, OptionState}
    composite_option = OptionState(
        option_name=string(option_names.Composite),
        description=(
            "Effective viscosity includes Peierls creep, dislocation creep and diffusion creep "
            *"flow laws depending on user input parameters."
            )
        )
    Dict{Int, OptionState}(
        option_ids[option_names.XCoordinateDependent] =>
            OptionState(
                option_name=string(option_names.XCoordinateDependent),
                description=(
                    "Viscosity is a function of lateral position of marker and pressure gradient."
                    ),
            ),
        option_ids[option_names.Isoviscous] =>
            OptionState(
                option_name=string(option_names.Isoviscous),
                description="Viscosity is constant",
            ),
        option_ids[option_names.Composite] => composite_option,
        # Maintain backward compatibility with old models
        option_ids[option_names.Composite2] => composite_option,
        option_ids[option_names.TemperatureDependentCouetteBenchmark] =>
            OptionState(
                option_name=string(option_names.TemperatureDependentCouetteBenchmark),
                description="Temperature-dependent viscosity for Couette flow benchmark.",
            ),
        option_ids[option_names.TemperatureDependentConvectionBenchmark] =>
            OptionState(
                option_name=string(option_names.TemperatureDependentConvectionBenchmark),
                description="Temperature-dependent viscosity for convection in a box benchmark.",
            )
    )
end

end # module 