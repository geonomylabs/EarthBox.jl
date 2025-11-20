module PlotUtils

import CairoMakie
import ...PlotParametersManager: PlotParameters

function get_figsize_pixels(parameters::PlotParameters)::Tuple{Int, Int}
    figsize = parameters.image.figsize
    dpi = parameters.image.figure_dpi
    figsize_pixels = (figsize[1] * dpi, figsize[2] * dpi)
    return figsize_pixels
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