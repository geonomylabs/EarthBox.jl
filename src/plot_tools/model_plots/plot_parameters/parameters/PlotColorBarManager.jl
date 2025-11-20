module PlotColorBarManager

import CairoMakie
import ...PlotDtypes: PlotDictType, AbstractPlotParameterGroup

Base.@kwdef mutable struct PlotColorBar <: AbstractPlotParameterGroup
    axis_fraction_for_color_bar::Float64 = 0.0 # obsolete
    axis_fraction_for_color_bar_gap::Float64 = 0.0 # obsolete
    color_map::Union{String, Symbol} = "None"
    minimum_value::Float64 = 0.0
    maximum_value::Float64 = 0.0
end

function PlotColorBar(plot_dict::PlotDictType)
    params = plot_dict["general_parameters"]
    return PlotColorBar(
        axis_fraction_for_color_bar=get(params, "axis_fraction_for_color_bar", 0.0), # obsolete
        axis_fraction_for_color_bar_gap=get(params, "axis_fraction_for_color_bar_gap", 0.0), # obsolete
        color_map=get(params, "color_map", "None"),
        minimum_value=get(params, "minimum_value", 0.0),
        maximum_value=get(params, "maximum_value", 0.0),
    )
end

function plot_colorbar!(
    fig::CairoMakie.Figure,
    limits::Tuple{Float64, Float64},
    color_map::Symbol;
    irow::Int=1,
    icol::Int=2,
    label::String="",
    vertical::Bool=true
)::Nothing
    @assert irow >= 1 && icol >= 1 "irow and icol must be greater than or equal to 1"
    CairoMakie.Colorbar(
        fig[irow, icol], 
        colormap=color_map, 
        limits=limits,
        label=label_pretty_filter(label),
#        vertical=false,
        #labelsize=25,
        #ticklabelsize=15
        )
    return nothing
end

function label_pretty_filter(label::String)::String
    if label == "C"
        return "Temperature (°C)"
    elseif label == "K"
        return "Temperature (K)"
    elseif label == "cm/yr"
        return "Velocity (cm/yr)"
    elseif label == "log10(Pas)"
        return "Viscosity (log₁₀ Pa⋅s)"
    elseif label == "log10(1/s)"
        return "Strain Rate (log₁₀ 1/s)"
    end
    return label
end

end # module