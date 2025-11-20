module Charts

import Plots

function plot_ncurves(chart_input::Dict{String, Any})::Nothing
    plot_dir = dirname(chart_input["plot_file_path"])
    check_output_directory(plot_dir)
  
    dpi = get(chart_input, "figure_dpi", 150)
    figsize = get(chart_input, "figsize", (5, 5))
    figsize_pixels = (figsize[1] * dpi, figsize[2] * dpi)
   
    axis_labels = get(chart_input, "axis_labels", ["X Axis", "Y Axis"])
    
    p = Plots.plot(
        title=get(chart_input, "title", "Plot Title"),
        xlabel=axis_labels[1],
        ylabel=axis_labels[2],
        legend=get(chart_input, "legend_location", :topleft),
        aspect_ratio=get(chart_input, "aspect_ratio", :auto),
        #fontfamily=fontfamily,
        titlefontsize=get(chart_input, "titlefontsize", 15),
        legendfontsize=get(chart_input, "legendfontsize", 12),
        guidefontsize=get(chart_input, "guidefontsize", 12),
        tickfontsize=get(chart_input, "tickfontsize", 10),
        annotationfontsize=get(chart_input, "annotationfontsize", 12),
        size=figsize_pixels,
        margin=10Plots.mm
    )
    set_plot_axes!(p, chart_input)
    plot_all_curves!(p, chart_input)
    
    if chart_input["iuse_inversion"] == 1
        Plots.yflip!(p)
    end
    
    xtext = chart_input["boxtext_info"][1]
    ytext = chart_input["boxtext_info"][2]
    text = chart_input["boxtext_info"][3]
    Plots.annotate!(p, xtext, ytext, text)
   
    Plots.savefig(p, chart_input["plot_file_path"])
    return nothing
end

function check_output_directory(output_path::String)::Nothing
    if !isdir(output_path)
        println("Output directory not found. Creating directory: $output_path")
        mkpath(output_path)
    end
    nothing
end

function set_plot_axes!(p::Plots.Plot, chart_input::Dict{String, Any})::Nothing
    xmin = chart_input["plot_dimensions_xy"][1]
    xmax = chart_input["plot_dimensions_xy"][2]
    ymin = chart_input["plot_dimensions_xy"][3]
    ymax = chart_input["plot_dimensions_xy"][4]
    xtick_size = chart_input["xtick_size"]
    ytick_size = chart_input["ytick_size"]
    Plots.xlims!(p, (xmin, xmax))
    Plots.ylims!(p, (ymin, ymax))
    Plots.xticks!(p, xmin:xtick_size:xmax)
    Plots.yticks!(p, ymin:ytick_size:ymax)
    nothing
end

function plot_all_curves!(
    p::Plots.Plot,
    chart_input::Dict{String, Any}
)::Nothing
    x_arrays = chart_input["data_xy"]["x_arrays"]
    y_arrays = chart_input["data_xy"]["y_arrays"]
    labels = chart_input["data_xy"]["labels"]
    line_styles = chart_input["data_xy"]["line_styles"]
    colors = chart_input["data_xy"]["colors"]
    line_widths = chart_input["data_xy"]["line_widths"]
    marker_sizes = chart_input["data_xy"]["marker_sizes"]
    marker_edge_colors = chart_input["data_xy"]["marker_edge_colors"]
    marker_edge_widths = chart_input["data_xy"]["marker_edge_widths"]
    fill_styles = chart_input["data_xy"]["fill_styles"]
    line_colors = chart_input["data_xy"]["line_colors"]

    ncurves = length(x_arrays)
    for i in 1:ncurves
        Plots.plot!(
            p,
            x_arrays[i],
            y_arrays[i],
            label=labels[i],
            linecolor=line_colors[i],
            linestyle=line_styles[i],
            linewidth=line_widths[i],
            markersize=marker_sizes[i],
            markercolor=colors[i],
            markerstrokecolor=marker_edge_colors[i],
            markerstrokewidth=marker_edge_widths[i],
            marker=fill_styles[i]
        )
    end
    nothing
end

function make_plot_name(
    plot_base_name::String,
    itime_step::Int,
    extension::String=".png"
)::String
    return plot_base_name * "_" * string(itime_step) * extension
end

end # module 