module FinalizePlots

import CairoMakie
import Printf
import ...PlotParametersManager: PlotParameters
import ..PlotUtils: get_figsize_pixels, get_root_list, create_all_roots

function finalize_plot!(
    fig::CairoMakie.Figure,
    axes::CairoMakie.Axis,
    parameters::PlotParameters,
    base_name::String;
    units::Union{String, Nothing}=nothing,
    extension::String=".png"
)::Nothing
    set_plot_title!(axes, parameters, base_name, units)
    save_figure(fig, parameters, base_name; extension=extension)
    return nothing
end

function set_plot_title!(
    axes::CairoMakie.Axis,
    parameters::PlotParameters,
    base_name::String,
    units::Union{String, Nothing}
)::Nothing
    axes.title = make_title_string(parameters, base_name, units)
    axes.titlesize = parameters.fonts.title_fontsize
    return nothing
end

function make_title_string(
    parameters::PlotParameters,
    base_name::String,
    units::Union{String, Nothing}
)::String
    plot_time = parameters.time.plot_time
    plot_time_units = parameters.time.plot_time_units
    if isnothing(units)
        units = ""
    end
    return Printf.@sprintf("%s %s : %10.4f %s", base_name, units, plot_time, plot_time_units)
end

function save_figure(
    fig::CairoMakie.Figure,
    parameters::PlotParameters, 
    base_name::String; 
    extension::String=".png"
)::Nothing
    plot_name = get_plot_name(parameters, base_name; extension=extension)
    
    # This is where most of the time is being spent!
    CairoMakie.save(plot_name, fig)

    return nothing
end

function get_plot_name(
    parameters::PlotParameters,
    base_name::String;
    extension::String=".png"
)::String
    ioutput = parameters.time.ioutput
    outpath = check_output_direc(parameters)
    plotname = joinpath(outpath, base_name * "_" * string(ioutput) * extension)
    return plotname
end

function check_output_direc(parameters::PlotParameters)::String
    outpath = parameters.paths.outpath
    outpath = check_direc(outpath)
    return outpath
end

function check_direc(outpath::String)::String
    if !isdir(outpath)
        println("Directory not found. Creating directory: $outpath")
        root_list = get_root_list(outpath)
        create_all_roots(root_list)
        mkdir(outpath)
    end
    return outpath
end

end