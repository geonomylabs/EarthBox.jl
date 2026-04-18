module PlotUtils

import CairoMakie
import ...PlotParametersManager: PlotParameters

function get_figsize_pixels(parameters::PlotParameters)::Tuple{Int, Int}
    figsize = parameters.image.figsize
    dpi = parameters.image.figure_dpi
    # Inch×dpi need not be integral; Cairo/Makie expect integer pixel sizes.
    w = round(Int, figsize[1] * dpi)
    h = round(Int, figsize[2] * dpi)
    return (w, h)
end

"""
    figsize_for_domain_extent(figsize_inches, dimensions; axis_width_fraction=0.85) -> (w, h)

Return figure size in inches so that, with `DataAspect()`, the **main axis** (not
the full figure width) matches the domain aspect `(xmax-xmin)/(ymax-ymin)`.

Grid plots use two columns `[axis | colorbar]`; `DataAspect()` applies to the axis
cell, which is only a **fraction** of the figure width. Height must therefore use
`figwidth * axis_width_fraction * dy/dx`, not `figwidth * dy/dx`, or the figure
row is too tall and Makie letterboxes (white bands above/below the heatmap).

The **figure width** `figsize_inches[1]` is kept; height is computed from that.
`axis_width_fraction` (default 0.85) approximates axis_column_width / figure_width;
tune if your layout differs (e.g. very wide colorbar → try 0.78).
"""
function figsize_for_domain_extent(
    figsize_inches::Tuple{Float64, Float64},
    dimensions::NTuple{4, Float64};
    axis_width_fraction::Float64 = 0.85,
)::Tuple{Float64, Float64}
    xmin, xmax, ymin, ymax = dimensions
    dx = xmax - xmin
    dy = ymax - ymin
    if dx <= 0 || dy <= 0
        return figsize_inches
    end
    w = figsize_inches[1]
    frac = clamp(axis_width_fraction, 0.5, 1.0)
    h = w * frac * dy / dx
    return (w, h)
end

function normalize_height_ratios(height_ratios::Vector{Float64})::Vector{Float64}
    total = sum(height_ratios)
    return [hr / total for hr in height_ratios]
end

function add_text_box(
    axes::CairoMakie.Axis,
    text::String,
    x::Float64,
    y::Float64;
    fontsize::Int=6
)::Nothing
    CairoMakie.text!(axes, x, y; text=text, align=(:left, :center), fontsize=fontsize)
    return nothing
end

function get_root_list(outpath::String)::Vector{String}
    outpath_root = dirname(outpath)
    outpath_root_root = dirname(outpath_root)
    outpath_root_root_root = dirname(outpath_root_root)
    root_list = [
        outpath_root_root_root,
        outpath_root_root,
        outpath_root
    ]
    return root_list
end

function create_all_roots(root_list::Vector{String})::Nothing
    for root in root_list
        if !isdir(root)
            println("Root directory not found. Creating directory: $root")
            mkpath(root)
        end
    end
    return nothing
end

end