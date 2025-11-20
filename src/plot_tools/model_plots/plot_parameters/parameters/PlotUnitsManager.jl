module PlotUnitsManager

import ...PlotDtypes: PlotDictType, AbstractPlotParameterGroup

Base.@kwdef mutable struct PlotUnits <: AbstractPlotParameterGroup
    active_units::Union{String, Nothing}
    length_units::String
    time_units::String
    velocity_units::String
    temperature_units::String
    viscosity_units::String
    stress_units::String
    strainrate_units::String
    pressure_units::String
end

function PlotUnits(plot_dict::PlotDictType)
    params = plot_dict["general_parameters"]
    return PlotUnits(
        active_units=nothing,
        length_units=get(params, "length_units", "None"),
        time_units=get(params, "time_units", "None"),
        velocity_units=get(params, "velocity_units", "None"),
        temperature_units=get(params, "temperature_units", "None"),
        viscosity_units=get(params, "viscosity_units", "None"),
        stress_units=get(params, "stress_units", "None"),
        strainrate_units=get(params, "strainrate_units", "None"),
        pressure_units=get(params, "pressure_units", "None"),
    )
end

function set_active_units!(units::PlotUnits, units_str::String)::Nothing
    units.active_units = units_str
    return nothing
end

end # module