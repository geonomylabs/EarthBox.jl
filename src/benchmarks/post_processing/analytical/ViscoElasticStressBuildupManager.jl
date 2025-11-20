module ViscoElasticStressBuildupManager

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
import ...Reader
import Printf: @sprintf

function compare_numerical_to_analytical(
    bench::Benchmarks
)::Tuple{Vector{Union{String, Float64}}, String}
    viscoelastic_stress_buildup = ViscoElasticStressBuildup(bench)
    make_benchmark_plot(viscoelastic_stress_buildup)
    return get_test_results(viscoelastic_stress_buildup)
end

mutable struct ViscoElasticStressBuildup
    itime_step::Int
    main_paths::Dict{String, String}
    test_info_dict::TestInfoDictType
    model_dict::ModelInputDictType
    materials_dict::Dict{Any, Any}
    materials_library_dict::MaterialsDictType
    relative_error_limit_percentage::Float64
    tkyr_vec::Vector{Float64}
    sxx_numerical::Vector{Float64}
    sxx_analytical::Vector{Float64}
end

function ViscoElasticStressBuildup(bench::Benchmarks)
    tkyr_vec, sxx_numerical, sxx_analytical = get_stress(bench)
    
    return ViscoElasticStressBuildup(
        bench.itime_step, bench.main_paths, bench.test_info_dict, bench.model_dict,
        bench.materials_dict, bench.materials_library_dict, 0.1,
        tkyr_vec, sxx_numerical, sxx_analytical
    )
end

function get_stress(
    bench::Benchmarks
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
    parameters_dict = make_parameters_dict(bench.model_dict)
    strain_rate = parameters_dict["strain_rate_bc"][1]

    nsteps, time_steps = get_time_steps(bench.itime_step)
    viscosity, shear_modulus = get_material_parameters(bench)

    tkyr_vec = zeros(Float64, nsteps)
    sxx_numerical = zeros(Float64, nsteps)
    sxx_analytical = zeros(Float64, nsteps)

    for itime_step in time_steps
        file_base_name = "Sxx_MPa"
        _gridx, _gridy, sxx, tmyr = Reader.read_scalar_grid(
            bench.main_paths["post_proc_input_path"], file_base_name, itime_step)
        sxx_numerical[itime_step] = maximum(sxx)
        tkyr_vec[itime_step] = tmyr * 1000.0
        tsec = tmyr * 1e6 * 365.0 * 24.0 * 3600.0
        sxx_analytical[itime_step] = calculate_analytical_stress(
            tsec, strain_rate, viscosity, shear_modulus
        )
    end
    return tkyr_vec, sxx_numerical, sxx_analytical
end

function get_material_parameters(
    bench::Benchmarks
)::Tuple{Float64, Float64}
    matid = 1
    material_name_layer1 = bench.materials_dict[matid]["mat_name"][1]
    material_dict_layer1 = bench.materials_library_dict[material_name_layer1]
    viscosity = material_dict_layer1["viscosity_iso"]
    shear_modulus = material_dict_layer1["shear_modulus"]
    return viscosity, shear_modulus
end

function get_time_steps(
    imax_time_step::Int
)::Tuple{Int, Vector{Int}}
    nsteps = imax_time_step
    time_steps = collect(1:nsteps)
    return nsteps, time_steps
end

function make_benchmark_plot(data::ViscoElasticStressBuildup)::Nothing
    data_xy = Dict{String, Any}(
        "x_arrays" => [
            data.tkyr_vec,
            data.tkyr_vec
        ],
        "y_arrays" => [
            data.sxx_analytical,
            data.sxx_numerical
        ],
        "labels" => ["Analytical", "Numerical"],
        "line_colors" => [:transparent, :transparent],
        "colors" => [:black, :red],
        "line_styles" => [:solid, :dash],
        "line_widths" => [1, 1],
        "marker_sizes" => [10, 5],
        "marker_edge_colors" => [:black, :black],
        "marker_edge_widths" => [0.5, 0.5],
        "fill_styles" => [:circle, :circle]
    )

    axis_labels = ["Time (Kyr)", "Normal Stress (MPa)"]
    plot_dimensions_xy = [0.0, 25.0, 0.0, 25.0]
    boxtext_x = plot_dimensions_xy[2] * 0.5
    boxtext_y = plot_dimensions_xy[4] * 0.05
    boxtext = ""
    boxtext_info = [boxtext_x, boxtext_y, boxtext]

    title_base_text = "Viscoelastic Stress Build-up Due to Pure Shear"
    title = "$(title_base_text)"

    plot_name = make_plot_name(
        data.main_paths["model_name"],
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
        "legend_location" => :bottomright,
        "legendfontsize" => 12,
        "figsize" => (10, 5),
        "guidefontsize" => 12,
        "titlefontsize" => 15,
        "tickfontsize" => 10,
        "annotation_fontsize" => 10,
        "xtick_size" => 5.0,
        "ytick_size" => 5.0
    )

    plot_ncurves(chart_input)
    return nothing
end

function get_test_results(
    data::ViscoElasticStressBuildup
)::Tuple{Vector{Union{String, Float64}}, String}
    result, result_msg = TestResults.get_test_results_numerical_vs_analytical(
        data.sxx_numerical,
        data.sxx_analytical,
        data.relative_error_limit_percentage / 100.0
    )
    return result, result_msg
end

""" Calculate analytical stress.
"""
function calculate_analytical_stress(
    time_seconds::Float64,
    strain_rate::Float64,
    viscosity::Float64,
    shear_modulus::Float64
)::Float64
    return (
        2.0 * strain_rate * viscosity
        * (1.0 - exp(-(time_seconds) * shear_modulus / viscosity)) / 1e6
    )
end

end # module 