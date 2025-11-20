module ChannelFlowVariableConductivityManager

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
    channel_flow_variable_conductivity = ChannelFlowVariableConductivity(bench)
    make_benchmark_plot(channel_flow_variable_conductivity)
    return get_test_results(channel_flow_variable_conductivity)
end

mutable struct ChannelFlowVariableConductivity
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

function ChannelFlowVariableConductivity(
    bench::Benchmarks
)
    xcoors_numerical_km, numerical_temperature_celsius, tmyr = 
        BenchmarkTools.get_xprofile_from_numerical_model(bench, "TempC")
    analytical_temperature_celsius = get_analytical_solution(
        xcoors_numerical_km, numerical_temperature_celsius, bench)
    
    return ChannelFlowVariableConductivity(
        bench.itime_step, bench.main_paths, bench.test_info_dict, bench.model_dict,
        bench.materials_dict, bench.materials_library_dict, tmyr, 2.2,
        xcoors_numerical_km, numerical_temperature_celsius, analytical_temperature_celsius
    )
end

function get_analytical_solution(
    xcoors_numerical_km::Vector{Float64},
    numerical_temperature_celsius::Vector{Float64},
    bench::Benchmarks
)::Vector{Float64}
    input_dict = get_input_for_analytical_calculation(
        xcoors_numerical_km, numerical_temperature_celsius, bench)
    analytical_temperature_kelvins = calculate_analytical_temperature(input_dict)
    return analytical_temperature_kelvins .- 273.0
end

function get_input_for_analytical_calculation(
    xcoors_numerical_km::Vector{Float64},
    numerical_temperature_celsius::Vector{Float64},
    bench::Benchmarks
)::Dict{String, Any}
    xnum = length(xcoors_numerical_km)
    width = xcoors_numerical_km[xnum] * 1000.0
    temperature_wall_kelvins = numerical_temperature_celsius[1] + 273.0

    parameters_dict = make_parameters_dict(bench.model_dict)
    ysize = parameters_dict["ysize"][1]
    ynum = parameters_dict["ynum"][1]
    dy_cell = ysize / (ynum - 1)
    height = ysize - dy_cell

    top_pressure_pa = parameters_dict["pressure_bc"][1]
    pressure_gradient = -top_pressure_pa / height

    material_name = bench.materials_dict[1]["mat_name"][1]
    material_dict = bench.materials_library_dict[material_name]
    viscosity = material_dict["viscosity_iso"]
    conductivity_ref = material_dict["thermal_conductivity_ref"]
    bcoef = material_dict["thermal_conductivity_a"]

    return Dict{String, Any}(
        "temperature_wall_kelvins" => temperature_wall_kelvins,
        "width" => width,
        "xnum" => xnum,
        "pressure_gradient" => pressure_gradient,
        "bcoef" => bcoef,
        "conductivity_ref" => conductivity_ref,
        "viscosity" => viscosity
    )
end

function make_benchmark_plot(data::ChannelFlowVariableConductivity)::Nothing
    data_xy = Dict{String, Any}(
        "x_arrays" => [
            data.xcoors_numerical_km,
            data.xcoors_numerical_km
        ],
        "y_arrays" => [
            data.analytical_temperature_celsius,
            data.numerical_temperature_celsius
        ],
        "labels" => ["Analytical Steady State", "Numerical Transient"],
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
    plot_dimensions_xy = [0.0, 30.0, 0.0, 900.0]
    boxtext_x = plot_dimensions_xy[2] * 0.75
    boxtext_y = plot_dimensions_xy[4] * 0.05
    boxtext = "Isoviscous\nShear Heating"
    boxtext_info = [boxtext_x, boxtext_y, boxtext]
    title_base_text = "Channel Flow with Variable Thermal Conductivity"
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
    data::ChannelFlowVariableConductivity
)::Tuple{Vector{Union{String, Float64}}, String}
    result, result_msg = TestResults.get_test_results_numerical_vs_analytical(
        data.numerical_temperature_celsius,
        data.analytical_temperature_celsius,
        data.relative_error_limit_percentage / 100.0
    )
    return result, result_msg
end

""" Calculate analytical temperature.
"""
function calculate_analytical_temperature(
    input_dict::Dict{String, Any}
)::Vector{Float64}
    temperature_wall_kelvins = input_dict["temperature_wall_kelvins"]
    xnum = input_dict["xnum"]
    viscosity = input_dict["viscosity"]
    conductivity_ref = input_dict["conductivity_ref"]
    bcoef = input_dict["bcoef"]
    pressure_gradient = input_dict["pressure_gradient"]
    width = input_dict["width"]
    
    temperature_kelvins = zeros(Float64, xnum)
    xmin = 0.0
    delta_x = width / (xnum - 1)
    for j in 1:xnum
        xcoor = xmin + delta_x * (j - 1)
        term1 = (
            width^4.0
            * bcoef / 192.0
            / conductivity_ref / temperature_wall_kelvins / viscosity
            * (pressure_gradient)^2.0
        )
        term2 = 1.0 - (2.0 * xcoor / width - 1.0)^4.0
        term3 = exp(term1 * term2)
        temperature_kelvins[j] = (
            temperature_wall_kelvins / bcoef * (term3 + bcoef - 1.0)
        )
    end

    return temperature_kelvins
end 

end # module
