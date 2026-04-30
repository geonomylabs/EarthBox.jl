module Charts

using CairoMakie

const LEGEND_POSITION_MAP = Dict{Symbol, Symbol}(
    :topleft     => :lt,
    :topright    => :rt,
    :bottomleft  => :lb,
    :bottomright => :rb,
    :bottom      => :cb,
    :top         => :ct,
    :left        => :lc,
    :right       => :rc,
)

function plot_ncurves(chart_input::Dict{String, Any})::Nothing
    plot_dir = dirname(chart_input["plot_file_path"])
    check_output_directory(plot_dir)

    dpi = get(chart_input, "figure_dpi", 150)
    figsize = get(chart_input, "figsize", (5, 5))
    figsize_pixels = (figsize[1] * dpi, figsize[2] * dpi)

    axis_labels = get(chart_input, "axis_labels", ["X Axis", "Y Axis"])

    fig = Figure(size = figsize_pixels)
    ax = Axis(
        fig[1, 1];
        title           = get(chart_input, "title", "Plot Title"),
        xlabel          = axis_labels[1],
        ylabel          = axis_labels[2],
        titlesize       = get(chart_input, "titlefontsize", 15),
        xlabelsize      = get(chart_input, "guidefontsize", 12),
        ylabelsize      = get(chart_input, "guidefontsize", 12),
        xticklabelsize  = get(chart_input, "tickfontsize", 10),
        yticklabelsize  = get(chart_input, "tickfontsize", 10),
    )

    aspect = _aspect_for(get(chart_input, "aspect_ratio", :auto))
    if aspect !== nothing
        ax.aspect = aspect
    end

    set_plot_axes!(ax, chart_input)
    plot_all_curves!(ax, chart_input)

    if chart_input["iuse_inversion"] == 1
        ax.yreversed = true
    end

    has_label = any(!isempty, chart_input["data_xy"]["labels"])
    if has_label
        position = get(LEGEND_POSITION_MAP, get(chart_input, "legend_location", :topleft), :lt)
        axislegend(ax; position = position,
                   labelsize = get(chart_input, "legendfontsize", 12))
    end

    xtext, ytext, text = chart_input["boxtext_info"]
    if !isempty(text)
        text!(ax, xtext, ytext;
              text = text,
              fontsize = get(chart_input, "annotationfontsize", 12))
    end

    save(chart_input["plot_file_path"], fig)
    return nothing
end

function check_output_directory(output_path::String)::Nothing
    if !isdir(output_path)
        println("Output directory not found. Creating directory: $output_path")
        mkpath(output_path)
    end
    nothing
end

function set_plot_axes!(ax::Axis, chart_input::Dict{String, Any})::Nothing
    xmin = chart_input["plot_dimensions_xy"][1]
    xmax = chart_input["plot_dimensions_xy"][2]
    ymin = chart_input["plot_dimensions_xy"][3]
    ymax = chart_input["plot_dimensions_xy"][4]
    xtick_size = chart_input["xtick_size"]
    ytick_size = chart_input["ytick_size"]
    xlims!(ax, xmin, xmax)
    ylims!(ax, ymin, ymax)
    ax.xticks = collect(xmin:xtick_size:xmax)
    ax.yticks = collect(ymin:ytick_size:ymax)
    nothing
end

function plot_all_curves!(ax::Axis, chart_input::Dict{String, Any})::Nothing
    d = chart_input["data_xy"]
    x_arrays            = d["x_arrays"]
    y_arrays            = d["y_arrays"]
    labels              = d["labels"]
    line_styles         = d["line_styles"]
    colors              = d["colors"]
    line_widths         = d["line_widths"]
    marker_sizes        = d["marker_sizes"]
    marker_edge_colors  = d["marker_edge_colors"]
    marker_edge_widths  = d["marker_edge_widths"]
    fill_styles         = d["fill_styles"]
    line_colors         = d["line_colors"]

    for i in eachindex(x_arrays)
        line_visible   = line_colors[i] !== :transparent
        marker_visible = fill_styles[i] !== :none
        label = isempty(labels[i]) ? nothing : labels[i]
        # Convention: transparent fill = open marker. The codebase often
        # leaves marker_edge_widths at 0.0 in that case (Plots rendered an
        # outline anyway); ensure the outline is visible under Makie too.
        stroke_w = (marker_visible && colors[i] === :transparent && marker_edge_widths[i] == 0.0) ?
                   1.0 : marker_edge_widths[i]

        if line_visible && marker_visible
            scatterlines!(
                ax, x_arrays[i], y_arrays[i];
                color       = line_colors[i],
                linestyle   = line_styles[i],
                linewidth   = line_widths[i],
                marker      = fill_styles[i],
                markercolor = colors[i],
                markersize  = marker_sizes[i],
                strokecolor = marker_edge_colors[i],
                strokewidth = stroke_w,
                label       = label,
            )
        elseif line_visible
            lines!(
                ax, x_arrays[i], y_arrays[i];
                color     = line_colors[i],
                linestyle = line_styles[i],
                linewidth = line_widths[i],
                label     = label,
            )
        elseif marker_visible
            scatter!(
                ax, x_arrays[i], y_arrays[i];
                marker      = fill_styles[i],
                color       = colors[i],
                markersize  = marker_sizes[i],
                strokecolor = marker_edge_colors[i],
                strokewidth = stroke_w,
                label       = label,
            )
        end
    end
    nothing
end

_aspect_for(value::Symbol) = value === :equal ? DataAspect() : nothing
_aspect_for(value::Real)   = AxisAspect(Float64(value))
_aspect_for(::Nothing)     = nothing

function make_plot_name(
    plot_base_name::String,
    itime_step::Int,
    extension::String=".png"
)::String
    return plot_base_name * "_" * string(itime_step) * extension
end

end # module
