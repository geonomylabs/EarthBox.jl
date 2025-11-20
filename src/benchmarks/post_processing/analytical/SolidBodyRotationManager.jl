module SolidBodyRotationManager

import EarthBox.InputTools.Reader: make_parameters_dict
import EarthBox.PlotToolsManager.Charts: plot_ncurves
import EarthBox.PlotToolsManager.Charts: make_plot_name
import EarthBox.EarthBoxDtypes: TestInfoDictType
import EarthBox.EarthBoxDtypes: ModelInputDictType
import EarthBox.EarthBoxDtypes: MaterialsDictType
import EarthBox.PlotSettingsManager: PLOT_SETTINGS
import ...BenchmarkTools
import ...TestResults
import ...BenchmarksStruct: Benchmarks
import ...Profiles: get_xprofile_scalar
import Printf: @sprintf

function compare_to_initial_temperature_wave(
    bench::Benchmarks
)::Tuple{Vector{Union{String, Float64}}, String}
    solid_body_rotation = SolidBodyRotation(bench)
    make_benchmark_plot(solid_body_rotation)
    return get_test_results(solid_body_rotation)
end

mutable struct SolidBodyRotation
    itime_step::Int
    main_paths::Dict{String, String}
    test_info_dict::TestInfoDictType
    model_dict::ModelInputDictType
    materials_dict::Dict{Any, Any}
    materials_library_dict::MaterialsDictType
    tmyr::Float64
    relative_error_limit_percentage::Float64
    xcoors_numerical_km::Vector{Float64}
    numerical_temperature_celsius_initial::Vector{Float64}
    numerical_temperature_celsius::Vector{Float64}
    subgrid_diff_coef_temp::Float64
    vy_cm_yr::Float64
    nrotations::Float64
end

function SolidBodyRotation(bench::Benchmarks)
    # Get initial temperature profile
    itime_step_initial = 1
    (
        _, numerical_temperature_celsius_initial, _
    ) = get_xprofile_scalar(
            bench.main_paths["post_proc_input_path"], 
            "TempC", itime_step_initial, 
            bench.test_info_dict[bench.main_paths["model_name"]]["y_index"]
            )
    
    # Get current temperature profile
    itime_step_current = bench.itime_step
    (
        xcoors_numerical_km, numerical_temperature_celsius, tmyr
    ) = 
        get_xprofile_scalar(
            bench.main_paths["post_proc_input_path"], 
            "TempC", itime_step_current, 
            bench.test_info_dict[bench.main_paths["model_name"]]["y_index"]
            )

    # Get parameters
    parameters_dict = make_parameters_dict(bench.model_dict)
    subgrid_diff_coef_temp = parameters_dict["subgrid_diff_coef_temp"][1]
    vy_cm_yr = parameters_dict["velocity_rotation"][1]
    
    nrotations = calculate_number_of_rotations(vy_cm_yr, xcoors_numerical_km, tmyr)
    println("nrotations: ", nrotations, " vy_cm_yr: ", vy_cm_yr, " tmyr: ", tmyr)
    return SolidBodyRotation(
        bench.itime_step, bench.main_paths, bench.test_info_dict, bench.model_dict,
        bench.materials_dict, bench.materials_library_dict, tmyr, 2.0,
        xcoors_numerical_km, numerical_temperature_celsius_initial, 
        numerical_temperature_celsius, subgrid_diff_coef_temp, vy_cm_yr, nrotations
    )
end

function calculate_number_of_rotations(
    vy_cm_yr::Float64,
    xcoors_numerical_km::Vector{Float64},
    tmyr::Float64
)::Float64
    vy_m_s = vy_cm_yr * 1e-2 / (365.0 * 24.0 * 60.0 * 60.0)
    xnum = length(xcoors_numerical_km)
    radius_m = xcoors_numerical_km[xnum] * 0.5 * 1000.0
    tseconds = tmyr * 1e6 * 365 * 24.0 * 60.0 * 60.0
    theta = vy_m_s / radius_m * tseconds
    nrotations = theta / (2.0 * π)
    return nrotations
end

function make_benchmark_plot(data::SolidBodyRotation)::Nothing
    data_xy = Dict{String, Any}(
        "x_arrays" => [
            data.xcoors_numerical_km,
            data.xcoors_numerical_km
        ],
        "y_arrays" => [
            data.numerical_temperature_celsius_initial,
            data.numerical_temperature_celsius
        ],
        "labels" => ["Numerical Initial", "Numerical Transient"],
        "line_colors" => [:blue, :red],
        "colors" => [:blue, :red],
        "line_styles" => [:dot, :dot],
        "line_colors" => [:blue, :transparent],
        "line_widths" => [1, 1],
        "marker_sizes" => [5, 5],
        "marker_edge_colors" => [:black, :black],
        "marker_edge_widths" => [0.5, 0.5],
        "fill_styles" => [:circle, :circle]
    )

    axis_labels = ["x (km)", "Temperature (Celsius)"]
    plot_dimensions_xy = [0.0, 50.0, 700.0, 1300.0]
    
    boxtext_x = plot_dimensions_xy[2] * 0.65
    boxtext_y = plot_dimensions_xy[4] * 0.98
    boxtext = @sprintf(
        "Number of Rotations: %.1f\nSubgrid Diff. Coefficient: %.2e",
        data.nrotations, data.subgrid_diff_coef_temp
    )
    boxtext_info = [boxtext_x, boxtext_y, boxtext]

    title_base_text = "Solid Body Rotation with Temperature Wave"
    title = "$(title_base_text): Numerical Time $(@sprintf("%.2f", data.tmyr)) Myr"

    plot_name = make_plot_name(
        data.main_paths["model_name"] * "_temperature",
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
        "legend_location" => :topleft,
        "legendfontsize" => 12,
        "figsize" => (8, 8),
        "guidefontsize" => 12,
        "titlefontsize" => 15,
        "tickfontsize" => 10,
        "annotation_fontsize" => 10,
        "xtick_size" => 5.0,
        "ytick_size" => 100.0
    )

    plot_ncurves(chart_input)
    return nothing
end

function get_test_results(
    data::SolidBodyRotation
)::Tuple{Vector{Union{String, Float64}}, String}
    result, result_msg = TestResults.get_test_results_numerical_vs_analytical(
        data.numerical_temperature_celsius,
        data.numerical_temperature_celsius_initial,
        data.relative_error_limit_percentage / 100.0
    )
    return result, result_msg
end

end # module 