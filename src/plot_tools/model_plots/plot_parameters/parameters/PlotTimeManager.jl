module PlotTimeManager

import ...PlotDtypes: PlotDictType, AbstractPlotParameterGroup

Base.@kwdef mutable struct PlotTime <: AbstractPlotParameterGroup
    ioutput::Int = 0
    tMyr::Float64 = 0.0
    plot_time::Float64 = 0.0
    plot_time_units::String = "None"
end

function set_plot_time_info!(plot_time::PlotTime, time::Float64, units::String)::Nothing
    plot_time.plot_time = time
    plot_time.plot_time_units = units
    return nothing
end

function PlotTime(plot_dict::PlotDictType)::PlotTime
    plot_params = plot_dict["general_parameters"]
    return PlotTime(
        ioutput = get(plot_params, "ioutput", 0),
        tMyr = get(plot_params, "tMyr", 0.0),
        plot_time = get(plot_params, "plot_time", 0.0),
        plot_time_units = get(plot_params, "plot_time_units", "None")
    )
end

end # module