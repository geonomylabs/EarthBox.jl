module ChannelFlowNonNewtonianManager

using Printf
import EarthBox.InputTools.Reader: make_parameters_dict
import EarthBox.PlotToolsManager.Charts: plot_ncurves
import EarthBox.PlotToolsManager.Charts: make_plot_name
import EarthBox.EarthBoxDtypes: MaterialsDictType
import EarthBox.EarthBoxDtypes: TestInfoDictType
import EarthBox.EarthBoxDtypes: ModelInputDictType
import EarthBox.PlotSettingsManager: PLOT_SETTINGS
import ...BenchmarkTools
import ...TestResults
import ...BenchmarksStruct: Benchmarks

function compare_numerical_to_analytical(
    bench::Benchmarks
)::Tuple{Vector{Union{String, Float64}}, String}
    channel_flow_non_newtonian = ChannelFlowNonNewtonian(bench)
    make_benchmark_velocity_plot(channel_flow_non_newtonian)
    make_benchmark_viscosity_plot(channel_flow_non_newtonian)
    return get_test_results(channel_flow_non_newtonian)
end

mutable struct ChannelFlowNonNewtonian
    itime_step::Int
    main_paths::Dict{String, String}
    test_info_dict::TestInfoDictType
    model_dict::ModelInputDictType
    materials_dict::Dict{Any, Any}
    materials_library_dict::MaterialsDictType
    tmyr::Float64
    relative_error_limit_percentage::Float64
    xcoors_velocity_numerical_km::Vector{Float64}
    numerical_velocity_cm_yr::Vector{Float64}
    xcoors_viscosity_numerical_km::Vector{Float64}
    numerical_viscosity_log_pas::Vector{Float64}
    analytical_velocity_cm_yr::Vector{Float64}
    analytical_viscosity_log_pas::Vector{Float64}
end

function ChannelFlowNonNewtonian(
    bench::Benchmarks
)
    xcoors_velocity_numerical_km, numerical_velocity_cm_yr, tmyr = 
        BenchmarkTools.get_xprofile_from_numerical_model(bench, "vel_cmyr")
    xcoors_viscosity_numerical_km, numerical_viscosity_log_pas, _ = 
        BenchmarkTools.get_xprofile_from_numerical_model(bench, "log_eta_Pas")
    analytical_velocity_cm_yr, analytical_viscosity_log_pas = 
        get_analytical_xprofiles(
            xcoors_velocity_numerical_km, xcoors_viscosity_numerical_km, bench)
    
    return ChannelFlowNonNewtonian(
        bench.itime_step, bench.main_paths, bench.test_info_dict, bench.model_dict,
        bench.materials_dict, bench.materials_library_dict, tmyr, 1.1,
        xcoors_velocity_numerical_km, numerical_velocity_cm_yr,
        xcoors_viscosity_numerical_km, numerical_viscosity_log_pas,
        analytical_velocity_cm_yr, analytical_viscosity_log_pas
    )
end

function get_analytical_xprofiles(
    xcoors_velocity_numerical_km::Vector{Float64},
    xcoors_viscosity_numerical_km::Vector{Float64},
    bench::Benchmarks
)::Tuple{Vector{Float64}, Vector{Float64}}
    input_dict = get_input_for_analytical_calculation(
        xcoors_velocity_numerical_km, xcoors_viscosity_numerical_km, bench)
    analytical_velocity_cm_yr, analytical_viscosity_log_pas = 
        calculate_analytical_velocity_and_viscosity_profiles(input_dict)
    return analytical_velocity_cm_yr, analytical_viscosity_log_pas
end

function get_input_for_analytical_calculation(
    xcoors_velocity_numerical_km::Vector{Float64},
    xcoors_viscosity_numerical_km::Vector{Float64},
    bench::Benchmarks
)::Dict{String, Any}
    xnum = length(xcoors_velocity_numerical_km)
    xdim = xcoors_viscosity_numerical_km[xnum] * 1000.0
    material_name = bench.materials_dict[1]["mat_name"][1]
    material_dict = bench.materials_library_dict[material_name]
    stress_exponent = material_dict["stress_exponent_n_dc"]
    pre_exponential_dc = material_dict["pre_exponential_dc"]
    pre_exponential_dc = pre_exponential_dc * 2.0 / (1e6^stress_exponent)
    
    parameters_dict = make_parameters_dict(bench.model_dict)
    top_pressure_pa = parameters_dict["pressure_bc"][1]
    ysize = parameters_dict["ysize"][1]
    ynum = parameters_dict["ynum"][1]
    dy_cell = ysize / (ynum - 1)
    ydim = ysize - dy_cell
    
    return Dict{String, Any}(
        "stress_exponent" => stress_exponent,
        "pre_exponential_dc" => pre_exponential_dc,
        "top_pressure_pa" => top_pressure_pa,
        "xnum" => xnum,
        "xdim" => xdim,
        "ydim" => ydim
    )
end

function make_benchmark_velocity_plot(data::ChannelFlowNonNewtonian)::Nothing
    data_xy = Dict{String, Any}(
        "x_arrays" => [
            data.xcoors_velocity_numerical_km,
            data.xcoors_velocity_numerical_km
        ],
        "y_arrays" => [
            data.analytical_velocity_cm_yr,
            data.numerical_velocity_cm_yr
        ],
        "labels" => ["Analytical", "Numerical"],
        "colors" => [:blue, :red],
        "line_styles" => [:solid, :dash],
        "line_widths" => [1, 1],
        "line_colors" => [:blue, :transparent],
        "marker_sizes" => [5, 5],
        "marker_edge_colors" => [:black, :black],
        "marker_edge_widths" => [0.5, 0.05],
        "fill_styles" => [:none, :circle]
    )

    axis_labels = ["x (km)", "Velocity_y (cm/yr)"]
    plot_dimensions_xy = [0.0, 10.0, 0.0, 60.0]
    boxtext_x = plot_dimensions_xy[2] * 0.5
    boxtext_y = plot_dimensions_xy[4] * 0.05
    boxtext = ""
    boxtext_info = [boxtext_x, boxtext_y, boxtext]
    title_base_text = "Non-Newtonian Channel Flow"
    title = "$(title_base_text): Numerical Time $(@sprintf("%.2f", data.tmyr)) Myr"
    plot_name = make_plot_name(
        data.main_paths["model_name"] * "_velocity",
        data.itime_step, PLOT_SETTINGS.plot_extension
    )
    plot_file_path = joinpath(
        data.main_paths["post_proc_output_path"], plot_name
    )
    
    chart_input = Dict{String, Any}(
        "plot_file_path" => plot_file_path,
        "title" => title,
        "axis_labels" => axis_labels,
        "plot_dimensions_xy" => plot_dimensions_xy,
        "data_xy" => data_xy,
        "boxtext_info" => boxtext_info,
        "iuse_inversion" => 0,
        "aspect_ratio" => :auto,
        "figure_dpi" => 150,
        "legend_location" => :bottom,
        "legendfontsize" => 12,
        "figsize" => (5, 5),
        "guidefontsize" => 12,
        "titlefontsize" => 15,
        "tickfontsize" => 10,
        "annotation_fontsize" => 10,
        "xtick_size" => 2.0,
        "ytick_size" => 5.0
    )
    
    plot_ncurves(chart_input)
    return nothing
end

function make_benchmark_viscosity_plot(data::ChannelFlowNonNewtonian)::Nothing
    data_xy = Dict{String, Any}(
        "x_arrays" => [
            data.xcoors_viscosity_numerical_km,
            data.xcoors_viscosity_numerical_km
        ],
        "y_arrays" => [
            data.analytical_viscosity_log_pas,
            data.numerical_viscosity_log_pas
        ],
        "labels" => ["Analytical", "Numerical"],
        "colors" => [:transparent, :red],
        "line_styles" => [:dot, :dot],
        "line_colors" => [:transparent, :transparent],
        "line_widths" => [1, 1],
        "marker_sizes" => [10, 5],
        "marker_edge_colors" => [:black, :black],
        "marker_edge_widths" => [0.0, 0.0],
        "fill_styles" => [:circle, :circle]
    )

    axis_labels = ["x (km)", "log10(Effective Viscosity) (Pa.s)"]
    plot_dimensions_xy = [0.0, 10.0, 18.0, 23.0]
    boxtext_x = plot_dimensions_xy[2] * 0.5
    boxtext_y = plot_dimensions_xy[4] * 0.05
    boxtext = ""
    boxtext_info = [boxtext_x, boxtext_y, boxtext]
    title_base_text = "Non-Newtonian Channel Flow"
    title = "$(title_base_text): Numerical Time $(@sprintf("%.2f", data.tmyr)) Myr"
    plot_name = make_plot_name(
        data.main_paths["model_name"] * "_viscosity",
        data.itime_step, PLOT_SETTINGS.plot_extension
    )
    plot_file_path = joinpath(
        data.main_paths["post_proc_output_path"], plot_name
    )
    
    chart_input = Dict{String, Any}(
        "plot_file_path" => plot_file_path,
        "title" => title,
        "axis_labels" => axis_labels,
        "plot_dimensions_xy" => plot_dimensions_xy,
        "data_xy" => data_xy,
        "boxtext_info" => boxtext_info,
        "iuse_inversion" => 0,
        "aspect_ratio" => :auto,
        "figure_dpi" => 150,
        "legend_location" => :bottom,
        "legendfontsize" => 12,
        "figsize" => (10, 5),
        "guidefontsize" => 12,
        "titlefontsize" => 15,
        "tickfontsize" => 10,
        "annotation_fontsize" => 10,
        "xtick_size" => 2.0,
        "ytick_size" => 0.5
    )
    
    plot_ncurves(chart_input)
    return nothing
end

function get_test_results(
    data::ChannelFlowNonNewtonian
)::Tuple{Vector{Union{String, Float64}}, String}
    result, result_msg = TestResults.get_test_results_numerical_vs_analytical(
        data.numerical_velocity_cm_yr,
        data.analytical_velocity_cm_yr,
        data.relative_error_limit_percentage / 100.0
    )
    return result, result_msg
end

function calculate_analytical_velocity_and_viscosity_profiles(
    input_dict::Dict{String, Any}
)::Tuple{Vector{Float64}, Vector{Float64}}
    stress_exponent = input_dict["stress_exponent"]
    pre_exponential_dc = input_dict["pre_exponential_dc"]
    top_pressure_pa = input_dict["top_pressure_pa"]
    xnum = input_dict["xnum"]
    xdim = input_dict["xdim"]
    ydim = input_dict["ydim"]

    cm_yr2m_s = 1.0 / (100.0 * 365.25 * 24.0 * 3600.0)
    pressure_gradient = -top_pressure_pa / ydim

    velocity_cm_yr = zeros(Float64, xnum)
    viscosity_log_pas = zeros(Float64, xnum)

    xmin = 0.0
    delta_x = xdim / (xnum - 1)
    for j in 1:xnum
        xcoor = xmin + delta_x * (j - 1)
        velocity = calc_velocity(
            pre_exponential_dc, stress_exponent, pressure_gradient,
            xdim, xcoor
        )
        velocity_cm_yr[j] = velocity / cm_yr2m_s
        if abs(xcoor - xdim / 2.0) > 0.0
            viscosity = calc_viscosity(
                pre_exponential_dc, stress_exponent, pressure_gradient,
                xdim, xcoor
            )
            viscosity_log_pas[j] = log10(viscosity)
        else
            viscosity_log_pas[j] = 0.0
        end
    end
    return velocity_cm_yr, viscosity_log_pas
end

function calc_velocity(
    pre_exponential_dc::Float64,
    stress_exponent::Float64,
    pressure_gradient::Float64,
    xdim::Float64,
    xcoor::Float64
)::Float64
    velocity = (
        pre_exponential_dc / (stress_exponent + 1.0)
        * (-pressure_gradient)^stress_exponent
        * (
            (xdim / 2.0)^(stress_exponent + 1.0)
            - (xcoor - xdim / 2.0)^(stress_exponent + 1.0)
        )
    )
    return velocity
end

function calc_viscosity(
    pre_exponential_dc::Float64,
    stress_exponent::Float64,
    pressure_gradient::Float64,
    xdim::Float64,
    xcoor::Float64
)::Float64
    viscosity = (
        1.0 / pre_exponential_dc
        * (-pressure_gradient)^(1.0 - stress_exponent)
        * (xcoor - xdim / 2.0)^(1.0 - stress_exponent)
    )
    return viscosity
end 

end # module