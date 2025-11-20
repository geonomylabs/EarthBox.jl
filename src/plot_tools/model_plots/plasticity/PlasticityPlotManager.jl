module PlasticityPlotManager

import EarthBox.ModelDataContainer: ModelData
import ...Charts: plot_ncurves, check_output_directory
import ..PlotDtypes: AxesType
import DataStructures: OrderedDict
using YAML
using Plots
using Printf

function get_stokes_convergence_plot_args_string()::String
    return """
# Stokes Convergence PlotKeyword Arguments
- `figsize::Tuple{Float64, Float64}`:
    - Figure size in inches in x-y directions. Default is (80.0, 10.0).
- `xspacing::Int`:
    - x-axis spacing. Default is 5.
- `log_l2_ymin::Int`:
    - Minimum y value for log L2 plot. Default is -5.
- `log_l2_ymax::Int`:
    - Maximum y value for log L2 plot. Default is 8.
- `log_l2_yspacing::Int`:
    - Y-axis spacing for log L2 plot. Default is 1.
- `plot_title::String`:
    - Plot title. Default is "Plasticity Iterations".
- `legend_fontsize::Int`:
    - Legend font size. Default is 8.
- `axis_label_size::Int`:
    - Axis label size. Default is 12.
- `axis_title_size::Int`:
    - Axis title size. Default is 14.
- `xtick_label_size::Int`:
    - X-tick label size. Default is 7. Currently not used.
- `ytick_label_size::Int`:
    - Y-tick label size. Default is 12.
- `annotation_font_size::Int`:
    - Annotation font size. Default is 4.
- `plot_yield_error::Bool`:
    - Whether to plot yield error. Default is false.
- `extension::String`:
    - Output file extension. Default is ".png".
"""
end

""" 
    run_plasticity_plotter(file_dir::String; kwargs...)

Run plasticity plotter.

# Arguments
- `file_dir::String`:
    - Directory containing the convergence.yml file.

$(get_stokes_convergence_plot_args_string())

"""
function run_plasticity_plotter(file_dir::String; kwargs...)
    file_name = "convergence.yml"
    file_path = joinpath(file_dir, file_name)
    data_dict = read_yaml_file(file_path)
    #print_data_dict(data_dict)
    curve_dict = initialize_curve_dict()
    fill_curve_dictionary(data_dict, curve_dict)
    make_plasticity_plot(file_dir, curve_dict; kwargs...)
end

""" Initialize curve dictionary for storing Picard iteration data.

Each key in the dictionary is a vector of values for each Picard iteration.
Keys refer to particular categories of data.
"""
function initialize_curve_dict()::Dict{String, Vector{Any}}
    curve_dict = Dict{String, Vector{Any}}()
    curve_dict["iter_ids"] = Int[]
    curve_dict["iglobals"] = Int[]
    curve_dict["l2_norms"] = Float64[]
    curve_dict["inf_norms"] = Float64[]
    curve_dict["residual_nonlinear_l2_norms"] = Float64[]
    curve_dict["pr_l2_norms"] = Float64[]
    curve_dict["solu_l2_norms"] = Float64[]
    curve_dict["times_myr"] = Float64[]
    curve_dict["status_all"] = String[]
    curve_dict["last_ids"] = Int[]
    curve_dict["last_l2_norms"] = Float64[]
    curve_dict["last_time_myr"] = Float64[]
    curve_dict["last_success_ids"] = Int[]
    curve_dict["last_failure_ids"] = Int[]
    curve_dict["last_success"] = Float64[]
    curve_dict["last_failure"] = Float64[]
    curve_dict["yield_error_log10"] = Float64[]
    return curve_dict
end

""" Fill curve dictionary using the data dictionary read from yaml output file.

# Arguments
- `data_dict::Dict{String, Vector{Any}}`: Dictionary with Picard iteration data. 
  The first three rows are header information.
- `curve_dict::Dict{String, Vector{Any}}`: Dictionary to fill with processed data
"""
function fill_curve_dictionary(
    data_dict::OrderedDict{Any, Any},
    curve_dict::Dict{String, Vector{Any}}
)
    (success_str, failure_str) = get_string_msg()

    keys_list = collect(keys(data_dict))
    nkeys = length(keys_list)
    for i in 1:nkeys
        if i > 3  # Skip first 3 rows since these are header rows
            key = keys_list[i]
            l2norm_velocity = data_dict[key][2]
            inf_norm_velocity = data_dict[key][3]
            residual_nonlinear = data_dict[key][4]
            l2norm_pr = data_dict[key][5]
            l2norm_solu = data_dict[key][6]
            yield_error = data_dict[key][7]
            time_myr = data_dict[key][8]
            status = data_dict[key][9]

            l2norm_velocity_log = calc_log10(l2norm_velocity)
            inf_norm_velocity_log = calc_log10(inf_norm_velocity)
            residual_nonlinear_log = calc_log10(residual_nonlinear)
            l2norm_pr_log = calc_log10(l2norm_pr)
            l2norm_solu_log = calc_log10(l2norm_solu)

            yield_error_log10 = calculate_log_yield_error(yield_error)

            push!(curve_dict["yield_error_log10"], yield_error_log10)
            push!(curve_dict["iter_ids"], key)
            push!(curve_dict["iglobals"], data_dict[key][1])
            push!(curve_dict["l2_norms"], l2norm_velocity_log)
            push!(curve_dict["inf_norms"], inf_norm_velocity_log)
            push!(curve_dict["residual_nonlinear_l2_norms"], 
                    residual_nonlinear_log)
            push!(curve_dict["pr_l2_norms"], l2norm_pr_log)
            push!(curve_dict["solu_l2_norms"], l2norm_solu_log)
            push!(curve_dict["times_myr"], data_dict[key][8])
            push!(curve_dict["status_all"], status)

            if status == success_str
                push!(curve_dict["last_success_ids"], key)
                push!(curve_dict["last_success"], l2norm_velocity_log)
            elseif status == failure_str
                push!(curve_dict["last_failure_ids"], key)
                push!(curve_dict["last_failure"], l2norm_velocity_log)
            end

            if status in [success_str, failure_str]
                push!(curve_dict["last_ids"], key)
                push!(curve_dict["last_l2_norms"], l2norm_velocity_log)
                push!(curve_dict["last_time_myr"], time_myr)
            end
        end
    end
end

function get_string_msg()::Tuple{String, String}
    success_str = "Convergence_Successful"
    failure_str = "Convergence_Failed_(iteration_limit_reached)"
    return success_str, failure_str
end

function calculate_log_yield_error(yield_error::Float64)::Float64
    yield_error_log10 = 0.0
    try
        yield_error_log10 = log10(yield_error)
    catch
        yield_error_log10 = 0.0
    end
    return yield_error_log10
end

""" Make plot of plasticity convergence criteria.

# Arguments
- `file_dir::String`: Output directory for the plot
- `curve_dict::Dict{String, Vector{Any}}`: Dictionary containing curve data

# Keyword Arguments
See `run_plasticity_plotter` for available options.
"""
function make_plasticity_plot(
    file_dir::String,
    curve_dict::Dict{String, Vector{Any}};
    kwargs...
)::Nothing
    figsize = get(kwargs, :figsize, (80, 10))
    xmin = get(kwargs, :xmin, 0)
    xmax = get(kwargs, :xmax, calc_xmax(curve_dict))
    xspacing = get(kwargs, :xspacing, 5)
    log_l2_ymin = get(kwargs, :log_l2_ymin, -5)
    log_l2_ymax = get(kwargs, :log_l2_ymax, 8)
    log_l2_yspacing = get(kwargs, :log_l2_yspacing, 1)
    plot_title = get(kwargs, :plot_title, "Plasticity Iterations")
    legend_font_size = get(kwargs, :legend_font_size, 8)
    axis_label_size = get(kwargs, :axis_label_size, 12)
    axis_title_size = get(kwargs, :axis_title_size, 14)
    xtick_label_size = get(kwargs, :xtick_label_size, 7) # Currently not used
    ytick_label_size = get(kwargs, :ytick_label_size, 12)
    annotation_font_size = get(kwargs, :annotation_font_size, 4)
    plot_yield_error = get(kwargs, :plot_yield_error, false)
    extension = get(kwargs, :extension, ".png")

    data_xy = make_data_dict(curve_dict, plot_yield_error)

    axis_labels = ["Iteration #", "Relative L2 Norm or Log(Yield Error)"]

    plot_name = "convergence" * extension
    plot_file_path = joinpath(file_dir, plot_name)

    chart_input = Dict{String, Any}(
        "plot_dimensions_xy" => [xmin, xmax, log_l2_ymin, log_l2_ymax],
        "plot_file_path" => plot_file_path,
        "title" => plot_title,
        "axis_labels" => axis_labels,
        "plot_dimensions_xy" => [xmin, xmax, log_l2_ymin, log_l2_ymax],
        "data_xy" => data_xy,
        "boxtext_info" => [0.0, 0.0, ""],
        "iuse_inversion" => 0,
        "aspect_ratio" => :auto,
        "figure_dpi" => 150,
        "legend_location" => :bottomright,
        "legendfontsize" => legend_font_size,
        "figsize" => figsize,
        "guidefontsize" => axis_label_size,
        "titlefontsize" => axis_title_size,
        "tickfontsize" => ytick_label_size,
        "annotation_fontsize" => annotation_font_size,
        "xtick_size" => xspacing,
        "ytick_size" => log_l2_yspacing
    )

    plot_ncurves(chart_input)

    annotate_final_sub_iteration_with_time(curve_dict, annotation_font_size)
    return nothing
end

function calc_xmax(curve_dict::Dict{String, Vector{Any}})::Int
    iter_ids = curve_dict["iter_ids"]
    return maximum(iter_ids)
end

function calc_log10(value::Float64)::Float64
    if value > 0
        value_log = log10(value)
    else
        value_log = 0.0
    end
    return value_log
end

function read_yaml_file(yaml_file_path::String)::OrderedDict{Any, Any}
    data_dict = YAML.load_file(yaml_file_path; dicttype=OrderedDict{Any,Any})
    return data_dict
end

function print_data_dict(data_dict::OrderedDict{Any, Any})::Nothing
    for (key, value) in data_dict
        println("Key: ", key)
        if isa(value, AbstractArray) && length(value) > 5
            # For long arrays, show first few elements
            println("Value (first 5 elements): ", value[1:5], " ...")
        else
            println("Value: ", value)
        end
    end
    return nothing
end

"""
    make_data_dict(curve_dict::Dict{String, Vector{Any}}, 
                   plot_yield_error::Bool) -> Dict{String, Vector{Any}}

Make data_xy dictionary for plotting function.

# Arguments
- `curve_dict::Dict{String, Vector{Any}}`: Dictionary containing curve data
- `plot_yield_error::Bool`: Whether to include yield error in plot

# Returns
- `Dict{String, Vector{Any}}`: Dictionary formatted for plotting
"""
function make_data_dict(
    curve_dict::Dict{String, Vector{Any}},
    plot_yield_error::Bool
)::Dict{String, Vector{Any}}
    iter_ids = curve_dict["iter_ids"]
    success_ids = curve_dict["last_success_ids"]
    success = curve_dict["last_success"]
    failure_ids = curve_dict["last_failure_ids"]
    failure = curve_dict["last_failure"]
    l2_norms = curve_dict["l2_norms"]
    l2_pr = curve_dict["pr_l2_norms"]
    l2_resnl = curve_dict["residual_nonlinear_l2_norms"]
    yield_error = curve_dict["yield_error_log10"]

    data_xy = Dict{String, Vector{Any}}(
        "x_arrays" => [
            iter_ids,
            iter_ids,
            success_ids,
            failure_ids,
            iter_ids,
            iter_ids
        ],
        "y_arrays" => [
            l2_norms,
            l2_norms,
            success,
            failure,
            l2_pr,
            l2_resnl
        ],
        "labels" => [
            "Velocity",
            "Iterations",
            "Success",
            "Failure",
            "Pressure",
            "Nonlinear Residual"
        ],
        "colors" => ["black", "black", "green", "red", "orange", "green"],
        "line_colors" => ["black", "black", :transparent, :transparent, "orange", "green"],
        "line_styles" => [:solid, :solid, :dot, :dot, :solid, :solid],
        "line_widths" => [2.0, 2.0, 1.0, 1.0, 2.0, 2.0],
        "marker_sizes" => [0.0, 0.0, 6.0, 6.0, 0.0, 0.0],
        "marker_edge_colors" => [:black, :black, :black, :black, :black, :black],
        "marker_edge_widths" => [4.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        "fill_styles" => [:circle, :circle, :circle, :circle, :circle, :circle]
    )

    if plot_yield_error
        push!(data_xy["x_arrays"], iter_ids)
        push!(data_xy["y_arrays"], yield_error)
        push!(data_xy["labels"], "log(Yield Error)")
        push!(data_xy["colors"], "blue")
        push!(data_xy["line_colors"], "blue")
        push!(data_xy["line_styles"], :solid)
        push!(data_xy["line_widths"], 1.0)
        push!(data_xy["marker_sizes"], 6.0)
        push!(data_xy["marker_styles"], :none)
        push!(data_xy["marker_edge_colors"], "black")
        push!(data_xy["marker_edge_widths"], 1.0)
        push!(data_xy["fill_styles"], "full")
    end

    return data_xy
end

"""
    annotate_final_sub_iteration_with_time(curve_dict::Dict{String, Vector{Any}}, 
                                          annotation_font_size::Int)

Annotate final sub-iteration with time in Myr.
"""
function annotate_final_sub_iteration_with_time(
    curve_dict::Dict{String, Vector{Any}},
    annotation_font_size::Int
)
    for i in 1:length(curve_dict["last_time_myr"])
        tmyr = curve_dict["last_time_myr"][i]
        lastid = curve_dict["last_ids"][i]
        l2norm = curve_dict["last_l2_norms"][i]
        label = string(tmyr)
        Plots.annotate!(
            (lastid, l2norm),
            text(label, annotation_font_size, :center, :bottom)
        )
    end
end

end  # module
