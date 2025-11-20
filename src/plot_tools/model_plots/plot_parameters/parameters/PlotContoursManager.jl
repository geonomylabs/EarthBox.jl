module PlotContoursManager

import ...PlotDtypes: PlotDictType, PlotParametersType, AbstractPlotParameterGroup
import CairoMakie
import Printf

const AxesType = CairoMakie.Axis

Base.@kwdef mutable struct PlotContours <: AbstractPlotParameterGroup
    iplot_contours::Int = 0
    contour_interval::Float64 = 0.0
    excluded_value::Float64 = -1e38
    excluded_vals::Vector{Float64} = [-1e38]
    V::Vector{Float64} = [0.0]
    V2::Vector{Float64} = [0.0]
    linewidths::Vector{Float64} = [1.0]
    contour_description::String = ""
    iplot_contour_labels::Int = 0
    xy_location_contour_legend::Tuple{Float64, Float64} = (0.0, 0.0)
end

function PlotContours(plot_dict::PlotDictType)::PlotContours
    plot_params = plot_dict["general_parameters"]
    return PlotContours(
        iplot_contours = 0,
        contour_interval = 0.0,
        excluded_value = -1e38,
        excluded_vals = [-1e38],
        V = [0.0],
        V2 = [0.0],
        linewidths = [1.0],
        contour_description = "",
        iplot_contour_labels = plot_params["iplot_contour_labels"],
        xy_location_contour_legend = plot_params["xy_location_contour_legend"]
    )
end

function reset_contour_description!(contours::PlotContours)::Nothing
    contours.contour_description = "Marker Contour Line Descriptions\n"
    return nothing
end

function add_to_contour_description!(contours::PlotContours, new_description::String)::Nothing
    suffix = get_suffix_for_description(contours)
    contours.contour_description *= new_description * suffix
    return nothing
end

function get_suffix_for_description(contours::PlotContours)::String
    if contours.contour_description === nothing
        suffix = ""
    else
        suffix = "\n"
    end
    return suffix
end

function activate_contours!(contours::PlotContours)::Nothing
    contours.iplot_contours = 1
    return nothing
end

function update_contour_levels!(
    contours::PlotContours,
    value_min::Float64,
    value_max::Float64
)::Nothing
    contours.V = get_values_for_contour_lines(
        value_min, value_max,
        contours.contour_interval, contours.excluded_vals
    )
    contours.V2 = get_values_for_contour_lines(
        value_min, value_max,
        contours.contour_interval * 2.0, [-1e38]
    )
    return nothing
end

function get_values_for_contour_lines(
    scalar_min::Float64,
    scalar_max::Float64,
    contour_interval::Float64,
    excluded_vals::Vector{Float64}
)::Vector{Float64}
    nlevels = count_contour_levels(scalar_min, scalar_max, contour_interval)
    contour_values = Float64[]
    for i in 1:nlevels
        value = get_contour_value_from_level(scalar_min, contour_interval, i-1)
        if !(value in excluded_vals)
            push!(contour_values, value)
        end
    end
    return contour_values
end

function update_linewidths!(
    contours::PlotContours,
    value_min::Float64,
    value_max::Float64,
    contour_line_width::Float64=0.5
)::Nothing
    contours.linewidths = get_linewidth_for_contour_lines(
        value_min, value_max, contours.contour_interval, contour_line_width, contours.excluded_vals)
    return nothing
end

function get_linewidth_for_contour_lines(
    scalar_min::Float64,
    scalar_max::Float64,
    contour_interval::Float64,
    contour_line_width::Float64=0.5,
    excluded_vals::Vector{Float64}=[-1e38]
)::Vector{Float64}
    nlevels = count_contour_levels(scalar_min, scalar_max, contour_interval)
    linewidths = Float64[]
    icount = 0
    for i in 1:nlevels
        val = get_contour_value_from_level(scalar_min, contour_line_width, i-1)
        if !(val in excluded_vals)
            if icount == 0
                lwv = contour_line_width
                icount += 1
            else
                lwv = contour_line_width + contour_line_width * 0.25
                icount = 0
            end
            push!(linewidths, lwv)
        end
    end
    return linewidths
end

function count_contour_levels(
    scalar_min::Float64,
    scalar_max::Float64,
    contour_interval::Float64
)::Int
    return floor(Int, (scalar_max - scalar_min)/contour_interval) + 1
end

function get_contour_value_from_level(
    scalar_min::Float64,
    contour_interval::Float64,
    ilevel::Int
)::Float64
    return scalar_min + contour_interval * ilevel
end

function plot_contours!(
    axes::CairoMakie.Axis,
    contours::PlotContours,
    gridx::Vector{Float64},
    gridy::Vector{Float64},
    grid_scalar::Matrix{Float64};
    color::Symbol=:black,
    labelsize::Int=5,
    number_format::String="%6.1f"
)::Nothing
    if contours.iplot_contours == 1
        feasible, message = check_contour_feasibility(grid_scalar, contours.V)
        if !feasible
            println("Warning: Skipping contours - $message")
            return nothing
        end
        # Note that with grid_scalar(i,j), i is in y-direction and j is in x-direction but
        # contour! expects gridx to be associated with the i-index and gridy to be associated
        # with the j-index. So a transpose is needed.
        min = minimum(grid_scalar)
        max = maximum(grid_scalar)
        for level in contours.V
            if level < min || level > max
                continue
            end
            # Get index of current level in contours.V
            level_idx = findfirst(x -> x == level, contours.V)
            # Set label_flag to false for every other contour
            label_flag = Bool(contours.iplot_contour_labels) && iseven(level_idx)
            ct = nothing
            try
                ct = CairoMakie.contour!(
                    axes,
                    gridx, 
                    gridy, 
                    grid_scalar',
                    levels=[level],
                    color=color,
                    linewidth=contours.linewidths[1],
                    labels=label_flag,
                    labelsize=labelsize,
                    labelformatter=create_formatter(number_format)
                    )
            catch e
                println("Error in contour! for level: $level")
                return nothing
            end
            ct.labelsize = labelsize
            #if label_flag
            #    label_contours!(axes, ct, level, color, number_format)
            #end
        end

    end
    return nothing
end

""" Very simple example of manual labeling of contours

"""
function label_contours!(ax, cplot, value, color, number_format::String="%6.1f")::Nothing
    fac = 0.4
    contour_points = cplot.plots[2][1][]
    println(">> contour_points: $contour_points")
    npts = length(contour_points)

    xmin = 1e39
    xmax = -1e39
    for (i, p) in enumerate(contour_points)
        xpt = p[1]
        if xpt < xmin
            xmin = xpt
        end
        if xpt > xmax
            xmax = xpt
        end
    end
    xfocus = xmin + fac * (xmax - xmin)
    println(">> xmin: $xmin, xmax: $xmax, xfocus: $xfocus")

    focus_index = -1
    dist_min = 1e39
    for (i, p) in enumerate(contour_points)
        xpt = p[1]
        dist = abs(xpt - xfocus)
        if dist < dist_min
            dist_min = dist
            focus_index = i
        end
    end
    p_focus = contour_points[focus_index]
    println(">> p_focus: $p_focus")
    CairoMakie.text!(ax, [(create_formatter(number_format)(value), p_focus)], align=(:center, :center), fontsize=10, color=color)
    return nothing
end

function create_formatter(format_string::String)::Function
    if !startswith(format_string, "%")
        format_string = "%" * format_string
    end
    # Pre-compile the format for better performance
    try
        fmt = Printf.Format(format_string)
        return function(value::Real)
            try
                return Printf.format(fmt, value)
            catch e
                return string(round(value, digits=2))
            end
        end
    catch e
        # If format string is invalid, return a simple formatter
        return function(value::Real)
            return string(round(value, digits=2))
        end
    end
end

function get_customized_contour_colors(
    colors::String
)::Union{Vector{Tuple{Float64, Float64, Float64}}, String}
    if colors == "cobalt_blue"
        colors = [(0.0, 0.2, 0.6)]
    elseif colors == "crimson_red"
        colors = [(0.6, 0.0, 0.0)]
    end
    return colors
end

function update_excluded_values_list!(contours::PlotContours)::Nothing
    contours.excluded_vals = [contours.excluded_value]
    return nothing
end

function check_contour_feasibility(
    grid_scalar::Matrix{Float64}, 
    levels::Vector{Float64}
)::Tuple{Bool, String}
    """
    Check if the scalar data has sufficient variation to draw meaningful contours.
    
    Returns:
        - feasible: Boolean indicating if contours can be drawn
        - message: Description of why contours are not feasible (if applicable)
    """
    scalar_min = minimum(grid_scalar)
    scalar_max = maximum(grid_scalar)
    scalar_range = scalar_max - scalar_min
    
    if scalar_range < 1e-10
        return false, "No variation in scalar data (range = $scalar_range)"
    end
    
    levels_in_range = filter(level -> scalar_min <= level <= scalar_max, levels)
    if isempty(levels_in_range)
        return false, "No contour levels fall within data range [$scalar_min, $scalar_max]"
    else
        level_spacing = abs(levels[2] - levels[1])
        if scalar_range < level_spacing * 0.1  # Require at least 10% of level spacing
            return false, "Data range ($scalar_range) too small relative to contour spacing ($level_spacing)"
        end
    end    
    return true, "Contours feasible"
end

end # module