module Options

import EarthBox.OptionTools: OptionState
import EarthBox.OptionTools: make_option_names

const option_ids = Dict{Symbol, Int}(
    :Liao14 => 0,
    :Uniform => 1,
    :TemperatureDependent => 2,
    :Naliboff17 => 3,
    :Adiabatic => 4,
    :SekiguchiWaples => 5
)

const option_names = make_option_names(option_ids)

function get_options()::Dict{Int, OptionState}
    Dict{Int, OptionState}(
        option_ids[option_names.Liao14] =>
            OptionState(
                option_name=string(option_names.Liao14),
                description="This thermal conductivity model type is a function of "
                *"temperature and pressure following the equation described "
                *"by Liao et al. (2014)."
            ),
        option_ids[option_names.Uniform] =>
            OptionState(
                option_name=string(option_names.Uniform),
                description="This thermal conductivity model type involves constant "
                *"thermal conductivity."
            ),
        option_ids[option_names.TemperatureDependent] =>
            OptionState(
                option_name=string(option_names.TemperatureDependent),
                description="This thermal conductivity model type is a function of "
                *"temperature and is used for the channel flow benchmark with "
                *"temperature-dependent thermal conductivity."
            ),
        option_ids[option_names.Naliboff17] =>
            OptionState(
                option_name=string(option_names.Naliboff17),
                description="Thermal conductivity model from Naliboff et al. (2017). "
                *"This thermal conductivity model uses constant values within material "
                *"domains and increases thermal conductivity to 39.25 W/m/K to mimic "
                *"a sub-lithospheric adiabatic gradient."
            ),
        option_ids[option_names.Adiabatic] =>
            OptionState(
                option_name=string(option_names.Adiabatic),
                description="Thermal conductivity model for simple lithosphere model. "
                *"This thermal conductivity model uses constant values within material "
                *"domains and increases thermal conductivity to mimic a sub-lithospheric "
                *"adiabatic gradient."
            ),
        option_ids[option_names.SekiguchiWaples] =>
            OptionState(
                option_name=string(option_names.SekiguchiWaples),
                description="Thermal conductivity model adapted from Sekiguchi (1984) "
                *"and modified by Doug Waples as summarized by Hantschel and Kauerauf "
                *"(2009). Thermal conductivity is a function of temperature and matrix "
                *"thermal conductivity at a reference temperature of 20 degrees Celsius."
            )
    )
end

end # module 