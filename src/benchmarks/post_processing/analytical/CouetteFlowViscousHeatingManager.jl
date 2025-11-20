module CouetteFlowViscousHeatingManager

import EarthBox.InputTools.Reader: make_parameters_dict
import EarthBox.PlotToolsManager.Charts: plot_ncurves
import EarthBox.PlotToolsManager.Charts: make_plot_name
import EarthBox.EarthBoxDtypes: TestInfoDictType
import EarthBox.EarthBoxDtypes: ModelInputDictType
import EarthBox.EarthBoxDtypes: MaterialsDictType
import EarthBox.MathTools: linear_interp_vals!
import EarthBox.PlotSettingsManager: PLOT_SETTINGS
import ...BenchmarkTools
import ...TestResults
import ...BenchmarksStruct: Benchmarks
import Printf: @sprintf

function compare_numerical_to_analytical(
    bench::Benchmarks
)::Tuple{Vector{Union{String, Float64}}, String}
    couette_flow_viscous_heating = CouetteFlowViscousHeating(bench)
    make_benchmark_plot(couette_flow_viscous_heating)
    return get_test_results(couette_flow_viscous_heating)
end

mutable struct CouetteFlowViscousHeating
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
    numerical_interpolated_temperature_celsius::Vector{Float64}
    xcoors_analytical_m::Vector{Float64}
    analytical_temperature_kelvins::Vector{Float64}
    brinkman::Float64
    phi1::Float64
end

function CouetteFlowViscousHeating(bench::Benchmarks)
    xcoors_numerical_km, numerical_temperature_celsius, tmyr = 
        BenchmarkTools.get_xprofile_from_numerical_model(bench, "TempC")
    
    (
        xcoors_analytical_m, 
        analytical_temperature_kelvins, 
        brinkman, 
        phi1
    ) = get_analytical_solution(
            bench, xcoors_numerical_km, numerical_temperature_celsius)
    
    numerical_interpolated_temperature_celsius = 
        interpolate_numerical_to_analytical_coors(
            xcoors_numerical_km, 
            numerical_temperature_celsius,
            xcoors_analytical_m
        )

    return CouetteFlowViscousHeating(
        bench.itime_step, bench.main_paths, bench.test_info_dict, bench.model_dict,
        bench.materials_dict, bench.materials_library_dict, tmyr, 0.3,
        xcoors_numerical_km, numerical_temperature_celsius, 
        numerical_interpolated_temperature_celsius,
        xcoors_analytical_m, analytical_temperature_kelvins,
        brinkman, phi1
    )
end

function get_analytical_solution(
    bench::Benchmarks,
    xcoors_numerical_km::Vector{Float64},
    numerical_temperature_celsius::Vector{Float64}
)::Tuple{Vector{Float64}, Vector{Float64}, Float64, Float64}
    input_dict = get_input_for_analytical_calculation(
        bench, xcoors_numerical_km, numerical_temperature_celsius
    )
    return calculate_analytical_solutions(input_dict)
end

function get_input_for_analytical_calculation(
    bench::Benchmarks,
    xcoors_numerical_km::Vector{Float64},
    numerical_temperature_celsius::Vector{Float64}
)::Dict{String, Any}
    xnum = length(xcoors_numerical_km)
    xnum = length(xcoors_numerical_km)
    temperature_right_wall_kelvins = numerical_temperature_celsius[xnum] + 273.0
    temperature_left_wall_kelvins = numerical_temperature_celsius[1] + 273.0
    width = xcoors_numerical_km[xnum] * 1000.0
    gas_constant = 8.314  # J/mol/K
    
    material_name = bench.materials_dict[1]["mat_name"][1]
    material_dict = bench.materials_library_dict[material_name]
    activation_energy = material_dict["activation_energy_td"] * 1000.0  # J/mol
    
    return Dict{String, Any}(
        "temperature_right_wall_kelvins" => temperature_right_wall_kelvins,
        "temperature_left_wall_kelvins" => temperature_left_wall_kelvins,
        "width" => width,
        "activation_energy" => activation_energy,
        "gas_constant" => gas_constant
    )
end

function interpolate_numerical_to_analytical_coors(
    xcoors_numerical_km::Vector{Float64},
    numerical_temperature_celsius::Vector{Float64},
    xcoors_analytical_m::Vector{Float64}
)::Vector{Float64}
    ncoors_analytical = length(xcoors_analytical_m)
    numerical_interpolated_temperature_kelvins = zeros(Float64, ncoors_analytical)
    linear_interp_vals!(
        xcoors_numerical_km,
        numerical_temperature_celsius,
        xcoors_analytical_m./1000.0,
        numerical_interpolated_temperature_kelvins
    )
    return numerical_interpolated_temperature_kelvins
end

function make_benchmark_plot(data::CouetteFlowViscousHeating)::Nothing
    data_xy = Dict{String, Any}(
        "x_arrays" => [
            data.xcoors_analytical_m ./ 1000.0,
            data.xcoors_numerical_km
        ],
        "y_arrays" => [
            data.analytical_temperature_kelvins .- 273.0,
            data.numerical_temperature_celsius
        ],
        "labels" => ["Analytical", "Numerical"],
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
    plot_dimensions_xy = [0.0, 30.0, 0.0, 1400.0]
    boxtext_x = plot_dimensions_xy[2] * 0.5
    boxtext_y = plot_dimensions_xy[4] * 0.05
    boxtext = @sprintf("Brinkman Number: %.3f\nTemperature-dependent Viscosity", data.brinkman)
    boxtext_info = [boxtext_x, boxtext_y, boxtext]
    
    title_base_text = "Couette Flow with Shear Heating"
    title = "$(title_base_text): Numerical Time $(@sprintf("%.2f", data.tmyr)) Myr"
    
    plot_name = make_plot_name(
        data.main_paths["model_name"],
        data.itime_step,
        PLOT_SETTINGS.plot_extension
    )
    plot_file_path = joinpath(
        data.main_paths["post_proc_output_path"],
        plot_name
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
        "legend_location" => :bottomleft,
        "legendfontsize" => 12,
        "figsize" => (10, 5),
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
    data::CouetteFlowViscousHeating
)::Tuple{Vector{Union{String, Float64}}, String}

    result, result_msg = TestResults.get_test_results_numerical_vs_analytical(
        data.numerical_interpolated_temperature_celsius .+ 273.0,
        data.analytical_temperature_kelvins,
        data.relative_error_limit_percentage / 100.0
    )
    return result, result_msg
end

function calculate_analytical_solutions(
    input_dict::Dict{String, Any}
)::Tuple{Vector{Float64}, Vector{Float64}, Float64, Float64}
    gas_constant = input_dict["gas_constant"]
    activation_energy = input_dict["activation_energy"]
    temperature_left_wall_kelvins = input_dict["temperature_left_wall_kelvins"]
    temperature_right_wall_kelvins = input_dict["temperature_right_wall_kelvins"]
    width = input_dict["width"]

    phi1 = calculate_phi1(
        activation_energy,
        gas_constant,
        temperature_left_wall_kelvins,
        temperature_right_wall_kelvins
    )

    bterm_min = 1e-6
    bterm_max = 10.0

    # Lower limit
    brinkman_min = calc_brinkman_number_from_bterm(bterm_min)
    phi1_min = calc_phi1_from_bterm_and_brinkman(bterm_min, brinkman_min)
    dphi1_min = phi1_min - phi1

    # Upper limit
    brinkman_max = calc_brinkman_number_from_bterm(bterm_max)
    phi1_max = calc_phi1_from_bterm_and_brinkman(bterm_max, brinkman_max)
    dphi1_max = phi1_max - phi1

    brinkman, phi1 = calc_brinkman_and_phi1_using_bisection(
        phi1, dphi1_min, dphi1_max, bterm_min, bterm_max
    )

    x_coors_m, temperature_kelvins = calculate_solutions_at_xcoors(
        phi1,
        brinkman,
        width,
        gas_constant,
        temperature_left_wall_kelvins,
        activation_energy
    )

    return x_coors_m, temperature_kelvins, brinkman, phi1
end

function calculate_phi1(
    activation_energy::Float64,
    gas_constant::Float64,
    temperature_left_wall_kelvins::Float64,
    temperature_right_wall_kelvins::Float64
)::Float64
    theta1 = (
        activation_energy *
        (temperature_right_wall_kelvins - temperature_left_wall_kelvins) /
        gas_constant /
        temperature_left_wall_kelvins^2
    )
    return exp(theta1)
end

function calc_brinkman_number_from_bterm(bterm::Float64)::Float64
    brinkman = (
        bterm^2.0 / 2.0 *
        (1.0 - ((exp(bterm) - 1.0) / (exp(bterm) + 1.0))^2.0)
    )
    return brinkman
end

function calc_phi1_from_bterm_and_brinkman(
    bterm::Float64,
    brinkman::Float64
)::Float64
    return bterm^2.0 / 2.0 / brinkman
end

function calc_brinkman_and_phi1_using_bisection(
    phi1::Float64,
    dphi1_min::Float64,
    dphi1_max::Float64,
    bterm_min::Float64,
    bterm_max::Float64
)::Tuple{Float64, Float64}
    brinkman_cur = 0.0
    phi1_cur = 0.0
    
    if dphi1_min <= 0 <= dphi1_max
        dphi1_cur = 1e6
        while abs(dphi1_cur) > 1e-6
            bterm_cur = (bterm_min + bterm_max) / 2
            brinkman_cur = calc_brinkman_number_from_bterm(bterm_cur)
            phi1_cur = calc_phi1_from_bterm_and_brinkman(bterm_cur, brinkman_cur)
            dphi1_cur = phi1_cur - phi1
            
            if dphi1_cur <= 0
                bterm_min = bterm_cur
            else
                bterm_max = bterm_cur
            end
        end
    end
    
    return brinkman_cur, phi1_cur
end

function calculate_solutions_at_xcoors(
    phi1::Float64,
    brinkman::Float64,
    width::Float64,
    gas_constant::Float64,
    temperature_left_wall_kelvins::Float64,
    activation_energy::Float64
)::Tuple{Vector{Float64}, Vector{Float64}}
    ncoors = 100
    x_coors_m = zeros(Float64, ncoors)
    temperature_kelvins = zeros(Float64, ncoors)
    delta_phi = (phi1 - 1.0) / (ncoors - 1)
    
    for j in 1:ncoors
        phi = 1.0 + (j - 1) * delta_phi
        theta = log(phi)
        x_coors_m[j] = calculate_x_coordinate(brinkman, phi1, phi, width)
        temperature_kelvins[j] = calculate_temperature(
            gas_constant,
            theta,
            temperature_left_wall_kelvins,
            activation_energy
        )
    end
    
    return x_coors_m, temperature_kelvins
end

function calculate_x_coordinate(
    brinkman::Float64,
    phi1::Float64,
    phi::Float64,
    width::Float64
)::Float64
    C = sqrt(2.0 * brinkman * (phi1 - phi))
    B = sqrt(2.0 * brinkman * phi1)
    D = sqrt(2.0 * brinkman * (phi1 - 1.0))
    term1 = 1.0 / sqrt(2.0 * brinkman * phi1)
    fac1 = (D + B) / (D - B)
    fac2 = (C - B) / (C + B)
    return width * term1 * log(fac1 * fac2)
end

function calculate_temperature(
    gas_constant::Float64,
    theta::Float64,
    temperature_left_wall_kelvins::Float64,
    activation_energy::Float64
)::Float64
    return (
        temperature_left_wall_kelvins +
        theta * gas_constant * temperature_left_wall_kelvins^2.0 /
        activation_energy
    )
end

end # module 