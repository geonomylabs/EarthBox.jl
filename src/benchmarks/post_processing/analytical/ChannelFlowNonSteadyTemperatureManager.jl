module ChannelFlowNonSteadyTemperatureManager

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
import Printf: @sprintf

function compare_numerical_to_analytical(
    bench::Benchmarks
)::Tuple{Vector{Union{String, Float64}}, String}
    channel_flow_non_steady_temperature = ChannelFlowNonSteadyTemperature(bench)
    make_benchmark_plot(channel_flow_non_steady_temperature)
    return get_test_results(channel_flow_non_steady_temperature)
end

mutable struct ChannelFlowNonSteadyTemperature
    itime_step::Int
    main_paths::Dict{String, String}
    test_info_dict::TestInfoDictType
    model_dict::ModelInputDictType
    materials_dict::Dict{Any, Any}
    materials_library_dict::MaterialsDictType
    tmyr::Float64
    relative_error_limit_percentage::Float64
    xcoors_numerical_km::Vector{Float64}
    numerical_temperature_celsius::Vector{Float64}
    analytical_temperature_celsius::Vector{Float64}
end

function ChannelFlowNonSteadyTemperature(
    bench::Benchmarks
)
    xcoors_numerical_km, numerical_temperature_celsius, tmyr = 
        BenchmarkTools.get_xprofile_from_numerical_model(bench, "TempC")
    analytical_temperature_celsius = get_analytical_solution(
        xcoors_numerical_km, numerical_temperature_celsius, tmyr, bench)
    
    return ChannelFlowNonSteadyTemperature(
        bench.itime_step, bench.main_paths, bench.test_info_dict, bench.model_dict,
        bench.materials_dict, bench.materials_library_dict, tmyr, 0.7,
        xcoors_numerical_km, numerical_temperature_celsius, analytical_temperature_celsius
    )
end

function get_analytical_solution(
    xcoors_numerical_km::Vector{Float64},
    numerical_temperature_celsius::Vector{Float64},
    tmyr::Float64,
    bench::Benchmarks
)::Vector{Float64}
    input_dict = get_input_for_analytical_calculation(
        xcoors_numerical_km, numerical_temperature_celsius, tmyr, bench)
    analytical_temperature_celsius = calculate_analytical_temperature(input_dict)
    return analytical_temperature_celsius
end

function get_input_for_analytical_calculation(
    xcoors_numerical_km::Vector{Float64},
    numerical_temperature_celsius::Vector{Float64},
    tmyr::Float64,
    bench::Benchmarks
)::Dict{String, Any}
    parameters_dict = make_parameters_dict(bench.model_dict)

    ysize = parameters_dict["ysize"][1]
    ynum = parameters_dict["ynum"][1]
    dy_cell = ysize / (ynum - 1)
    channel_section_length = ysize - dy_cell

    top_pressure_pa = parameters_dict["pressure_bc"][1]
    pressure_gradient = -top_pressure_pa / channel_section_length

    xnum = length(xcoors_numerical_km)
    width = xcoors_numerical_km[xnum] * 1000.0

    temperature_top_celsius = parameters_dict["temperature_top"][1]
    temperature_top_kelvins = temperature_top_celsius + 273.0

    temperature_bottom_celsius = parameters_dict["temperature_bottom"][1]
    temperature_bottom_kelvins = temperature_bottom_celsius + 273.0

    thermal_gradient_avg = (
        (temperature_bottom_kelvins - temperature_top_kelvins) / ysize
    )

    y_index = bench.test_info_dict[bench.main_paths["model_name"]]["y_index"]
    y_test = (y_index-1) * ysize / (ynum - 1)
    temperature_initial = temperature_top_kelvins + y_test * thermal_gradient_avg

    material_name = bench.materials_dict[1]["mat_name"][1]
    material_dict = bench.materials_library_dict[material_name]

    viscosity_pas = material_dict["viscosity_iso"]
    conductivity_w_m_k = material_dict["thermal_conductivity_ref"]
    density_kg_m3 = material_dict["standard_density"]
    heat_capacity_j_k_kg = material_dict["heat_capacity"]

    return Dict{String, Any}(
        "tmyr" => tmyr,
        "temperature_initial" => temperature_initial,
        "xnum" => xnum,
        "viscosity_pas" => viscosity_pas,
        "density_kg_m3" => density_kg_m3,
        "heat_capacity_j_k_kg" => heat_capacity_j_k_kg,
        "conductivity_w_m_k" => conductivity_w_m_k,
        "thermal_gradient_avg" => thermal_gradient_avg,
        "pressure_gradient" => pressure_gradient,
        "width" => width
    )
end

function make_benchmark_plot(data::ChannelFlowNonSteadyTemperature)::Nothing
    data_xy = Dict{String, Any}(
        "x_arrays" => [
            data.xcoors_numerical_km,
            data.xcoors_numerical_km
        ],
        "y_arrays" => [
            data.analytical_temperature_celsius,
            data.numerical_temperature_celsius
        ],
        "labels" => ["Analytical Transient", "Numerical Transient"],
        "line_colors" => [:blue, :transparent],
        "colors" => [:blue, :red],
        "line_styles" => [:solid, :dash],
        "line_widths" => [1, 1],
        "marker_sizes" => [5, 5],
        "marker_edge_colors" => [:black, :black],
        "marker_edge_widths" => [0.5, 0.5],
        "fill_styles" => [:none, :circle]
    )

    axis_labels = ["x (km)", "Temperature (Celsius)"]
    plot_dimensions_xy = [0.0, 30.0, 500.0, 1600.0]
    boxtext_x = plot_dimensions_xy[2] * 0.75
    boxtext_y = plot_dimensions_xy[4] * 0.05
    boxtext = "Isoviscous\n"
    boxtext_info = [boxtext_x, boxtext_y, boxtext]
    title_base_text = "Isoviscous Channel Flow with Non-Steady Temperature"
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
        "figsize" => (8, 5),
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
    data::ChannelFlowNonSteadyTemperature
)::Tuple{Vector{Union{String, Float64}}, String}
    result, result_msg = TestResults.get_test_results_numerical_vs_analytical(
        data.numerical_temperature_celsius,
        data.analytical_temperature_celsius,
        data.relative_error_limit_percentage / 100.0
    )
    return result, result_msg
end

""" Calculate analytical temperature for isoviscous channel flow.
"""
function calculate_analytical_temperature(
    input_dict::Dict{String, Any}
)::Vector{Float64}
    tmyr = input_dict["tmyr"]
    temperature_initial = input_dict["temperature_initial"]
    xnum = input_dict["xnum"]
    viscosity_pas = input_dict["viscosity_pas"]
    density_kg_m3 = input_dict["density_kg_m3"]
    heat_capacity_j_k_kg = input_dict["heat_capacity_j_k_kg"]
    conductivity_w_m_k = input_dict["conductivity_w_m_k"]
    thermal_gradient_avg = input_dict["thermal_gradient_avg"]
    pressure_gradient = input_dict["pressure_gradient"]
    width = input_dict["width"]

    tsec = tmyr * 1e6 * 365.0 * 24.0 * 60.0 * 60.0
    # added 1/2 factor to match figure 16.5 in Gerya 2010
    ep_term = -1.0 / 4.0 / viscosity_pas * pressure_gradient * thermal_gradient_avg
    kappa = conductivity_w_m_k / density_kg_m3 / heat_capacity_j_k_kg
    mmax = 100
    temperature_celsius = zeros(Float64, xnum)
    xmin = 0.0
    delta_x = width / (xnum - 1)

    for j in 1:xnum
        x_coor = xmin + delta_x * (j - 1)
        delta_temp = 0.0
        for m in 1:mmax
            mfloat = Float64(m-1)
            fac1 = (π * (2.0 * mfloat - 1))^2.0
            fac2 = (π * (2.0 * mfloat - 1))^3.0
            emt_term = (
                width * width
                * (1.0 - exp(-kappa * tsec / width / width * fac1)) / (kappa * fac1)
            )
            fm_term = -8.0 * ep_term * width * width / fac2
            delta_temp = (
                delta_temp
                + fm_term * emt_term
                * sin(π * (2.0 * mfloat - 1) * x_coor / width)
            )
        end
        temperature_celsius[j] = (temperature_initial + delta_temp) - 273.0
    end
    return temperature_celsius
end

end # module 