module RayleighTaylorInstabilityManager

import Printf: @sprintf
import DataStructures: OrderedDict
import EarthBox.InputTools.Reader: make_parameters_dict
import EarthBox.PlotToolsManager.Charts: plot_ncurves
import EarthBox.PlotToolsManager.Charts: make_plot_name
import EarthBox.EarthBoxDtypes: TestInfoDictType
import EarthBox.EarthBoxDtypes: ModelInputDictType
import EarthBox.EarthBoxDtypes: MaterialsDictType
import EarthBox.StaggeredGrid.Spacing: initialize_spacing_arrays
import EarthBox.StaggeredGrid.Spacing: update_basic_1Dgrid_spacing!
import EarthBox.Interpolation.GridToMarker: get_marker_value_from_basic_grid_without_mapping_input
import EarthBox.ConversionFuncs: get_factor_cm_yr_to_m_s
import ...BenchmarkTools
import ...TestResults
import ...BenchmarksStruct: Benchmarks
import ...Reader: read_vector_grid

function compare_numerical_to_analytical(
    bench::Benchmarks
)::Tuple{Vector{Union{String, Float64}}, String}
    rayleigh_taylor_instability = RayleighTaylorInstability(bench)
    return get_test_results(rayleigh_taylor_instability)
end

mutable struct RayleighTaylorInstability
    itime_step::Int
    main_paths::Dict{String, String}
    test_info_dict::TestInfoDictType
    model_dict::ModelInputDictType
    materials_dict::Dict{Any, Any}
    materials_library_dict::MaterialsDictType
    relative_error_limit_percentage::Float64
    vy_wave_numerical_array::Vector{Float64}
    vy_wave_analytical_array::Vector{Float64}
end

function RayleighTaylorInstability(bench::Benchmarks)
    relative_error_limit_percentage = 3.5
    vy_wave_numerical_array = zeros(Float64, 1)
    vy_wave_analytical_array = zeros(Float64, 1)

    (
        _growth_factor_k_numerical,
        vy_wave_numerical_array[1],
        _growth_factor_k_analytical,
        vy_wave_analytical_array[1]
    ) = rayleigh_taylor_instability_quantities(bench)

    return RayleighTaylorInstability(
        bench.itime_step, bench.main_paths, bench.test_info_dict, bench.model_dict,
        bench.materials_dict, bench.materials_library_dict,
        relative_error_limit_percentage, vy_wave_numerical_array, vy_wave_analytical_array
    )
end

function get_test_results(
    data::RayleighTaylorInstability
)::Tuple{Vector{Union{String, Float64}}, String}
    result, result_msg = TestResults.get_test_results_numerical_vs_analytical(
        data.vy_wave_numerical_array,
        data.vy_wave_analytical_array,
        data.relative_error_limit_percentage / 100.0
    )
    return result, result_msg
end

function rayleigh_taylor_instability_quantities(
    bench::Benchmarks
)::Tuple{Float64, Float64, Float64, Float64}
    """ Calculate quantities for benchmark using numerical results.

    This benchmark involves a Rayleigh-Taylor instability as described by
    Ramberg (1968). Quantities are calculated using numerical solutions for
    comparison with analytical values.
    """
    file_base_name = "vel_cmyr"
    (
        gridx_b_km,
        gridy_b_km,
        _velocity_x_cm_yr,
        velocity_y_cm_yr,
        _tmyr
    ) = read_vector_grid(
        bench.main_paths["post_proc_input_path"],
        file_base_name,
        bench.itime_step
    )

    cm_yr_to_m_s = get_factor_cm_yr_to_m_s()
    velocity_y = velocity_y_cm_yr .* cm_yr_to_m_s

    gridx_b = gridx_b_km .* 1000.0
    gridy_b = gridy_b_km .* 1000.0
    xnum = length(gridx_b_km)
    ynum = length(gridy_b_km)
    xsize = gridx_b[xnum]
    ysize = gridy_b[ynum]

    (
        density_layer1,
        density_layer2,
        viscosity_layer1,
        viscosity_layer2
    ) = get_material_parameters(bench.materials_dict, bench.materials_library_dict)

    (
        layer1_thickness,
        wavelength,
        amplitude_initial,
        gravity_y
    ) = get_model_input_parameters(bench.model_dict)

    layer2_thickness = ysize - layer1_thickness

    xwave, ywave = calculate_wave_peak_location(
        xsize, layer1_thickness, amplitude_initial
    )

    velocity_y_wave_numerical = calculate_velocity_at_central_peak_of_instability_wave(
        gridx_b, gridy_b, velocity_y, xwave, ywave
    )

    q_factor = (
        (density_layer1 - density_layer2)
        * layer2_thickness * gravity_y / 2.0 / viscosity_layer2
    )

    growth_factor_k_numerical = calculate_growth_factor(
        velocity_y_wave_numerical, amplitude_initial, q_factor
    )

    (
        growth_factor_k_analytical,
        velocity_y_wave_analytical
    ) = calculate_analytical_wave_velocity_and_growth_rate(
        layer1_thickness,
        layer2_thickness,
        density_layer1,
        density_layer2,
        viscosity_layer1,
        viscosity_layer2,
        amplitude_initial,
        wavelength
    )

    m_s_to_cm_yr = 1.0 / cm_yr_to_m_s
    b1, b2 = get_scaling_coefficients_rayleigh_taylor()
    output_dict = OrderedDict{String, Any}(
        "layer1_thickness" => layer1_thickness,
        "layer2_thickness" => layer2_thickness,
        "density_layer1" => density_layer1,
        "density_layer2" => density_layer2,
        "viscosity_layer1" => viscosity_layer1,
        "viscosity_layer2" => viscosity_layer2,
        "amplitude_initial" => amplitude_initial,
        "wavelength" => wavelength,
        "growth_factor_k_numerical" => growth_factor_k_numerical,
        "velocity_y_wave_numerical(cm/yr)" => velocity_y_wave_numerical * m_s_to_cm_yr,
        "growth_factor_k_analytical" => growth_factor_k_analytical,
        "velocity_y_wave_analytical(cm/yr)" => velocity_y_wave_analytical * m_s_to_cm_yr,
        "theta1" => 2.0 * π * layer1_thickness / wavelength,
        "b1" => b1,
        "b2" => b2,
        "b1*K+b2 (numerical)" => b1 * growth_factor_k_numerical + b2,
        "b1*K+b2 (analytical)" => b1 * growth_factor_k_analytical + b2
    )

    print_benchmark_results(output_dict)

    return (
        growth_factor_k_numerical,
        velocity_y_wave_numerical,
        growth_factor_k_analytical,
        velocity_y_wave_analytical
    )
end

function get_scaling_coefficients_rayleigh_taylor()::Tuple{Float64, Float64}
    """ Get scaling coefficients for benchmark plot.
    """
    b1 = 0.5
    b2 = 0.2
    return b1, b2
end

function get_material_parameters(
    materials_dict::Dict{Any, Any},
    materials_library_dict::MaterialsDictType
)::Tuple{Float64, Float64, Float64, Float64}
    """ Get material parameters from input files.
    """
    matid = 2
    material_name_layer1 = materials_dict[matid]["mat_name"][1]
    material_dict_layer1 = materials_library_dict[material_name_layer1]
    density_layer1 = material_dict_layer1["standard_density"]
    viscosity_layer1 = material_dict_layer1["viscosity_iso"]

    matid = 3
    material_name_layer2 = materials_dict[matid]["mat_name"][1]
    material_dict_layer2 = materials_library_dict[material_name_layer2]
    density_layer2 = material_dict_layer2["standard_density"]
    viscosity_layer2 = material_dict_layer2["viscosity_iso"]
    return density_layer1, density_layer2, viscosity_layer1, viscosity_layer2
end

function get_model_input_parameters(
    model_dict::ModelInputDictType
)::Tuple{Float64, Float64, Float64, Float64}
    """ Get input from model input file.
    """
    parameters_dict = make_parameters_dict(model_dict)
    layer1_thickness = parameters_dict["depth_interface_h1"][1]
    wave_length_lambda = parameters_dict["wave_length_lambda"][1]
    amplitude_initial = parameters_dict["amplitude_initial"][1]
    gravity_y = parameters_dict["gravity_y"][1]
    return layer1_thickness, wave_length_lambda, amplitude_initial, gravity_y
end

function calculate_velocity_at_central_peak_of_instability_wave(
    gridx_b::Vector{Float64},
    gridy_b::Vector{Float64},
    velocity_y::Matrix{Float64},
    xwave::Float64,
    ywave::Float64
)::Float64
    """ Calculate numerical velocity of wave instability.
    """
    xnum = length(gridx_b)
    ynum = length(gridy_b)

    xstp_b, ystp_b = initialize_spacing_arrays(gridx_b, gridy_b)
    update_basic_1Dgrid_spacing!(xnum, gridx_b, xstp_b)
    update_basic_1Dgrid_spacing!(ynum, gridy_b, ystp_b)

    velocity_y_wave = get_marker_value_from_basic_grid_without_mapping_input(
        gridx_b, gridy_b, xstp_b, ystp_b, velocity_y, xwave, ywave
    )

    return velocity_y_wave
end

function calculate_wave_peak_location(
    xsize::Float64,
    layer1_thickness::Float64,
    amplitude_initial::Float64
)::Tuple{Float64, Float64}
    """ Calculate x and y-location of initial peak of instability.
    """
    xwave = xsize / 2.0
    ywave = layer1_thickness - amplitude_initial
    return xwave, ywave
end

function calculate_growth_factor(
    velocity_y_wave::Float64,
    amplitude_initial::Float64,
    q_factor::Float64
)::Float64
    """ Calculate non-dimensional growth factor.
    """
    growth_factor_k = abs(velocity_y_wave) / amplitude_initial / q_factor
    return growth_factor_k
end

function calculate_analytical_wave_velocity_and_growth_rate(
    layer1_thickness::Float64,
    layer2_thickness::Float64,
    density_layer1::Float64,
    density_layer2::Float64,
    viscosity_layer1::Float64,
    viscosity_layer2::Float64,
    amplitude_initial::Float64,
    wavelength::Float64
)::Tuple{Float64, Float64}
    """ Calculate analytical velocity and growth rate for RT benchmark.
    """
    theta1 = 2.0 * π * layer1_thickness / wavelength
    theta2 = 2.0 * π * layer2_thickness / wavelength

    c11_term = calculate_c11_term(
        viscosity_layer1, viscosity_layer2, theta1, theta2
    )

    d12_term = calculate_d12_term(
        viscosity_layer1, viscosity_layer2, theta1, theta2
    )

    j21_term = calculate_j21_term(
        viscosity_layer1, viscosity_layer2, theta1, theta2
    )

    j22_term = calculate_j22_term(
        viscosity_layer1, viscosity_layer2, theta1, theta2
    )

    growth_factor_k = -d12_term / (c11_term * j22_term - d12_term * j21_term)

    velocity_y_wave = (
        -growth_factor_k * (density_layer1 - density_layer2)
        * layer2_thickness * 9.81 / 2.0 / viscosity_layer2 * amplitude_initial
    )

    return growth_factor_k, velocity_y_wave
end

function calculate_c11_term(
    viscosity_layer1::Float64,
    viscosity_layer2::Float64,
    theta1::Float64,
    theta2::Float64
)::Float64
    """ Calculate C11 term for Rayleigh-Taylor instability benchmark.
    """
    term1 = (
        viscosity_layer1 * 2.0 * theta1^2.0 / viscosity_layer2
        / (cosh(2.0 * theta1) - 1.0 - 2.0 * theta1^2.0)
    )
    term2 = (
        2 * theta2^2.0 / (cosh(2.0 * theta2) - 1.0 - 2.0 * theta2^2.0)
    )
    c11_term = term1 - term2
    return c11_term
end

function calculate_d12_term(
    viscosity_layer1::Float64,
    viscosity_layer2::Float64,
    theta1::Float64,
    theta2::Float64
)::Float64
    """ Calculate d12 term for Rayleigh-Taylor instability benchmark.
    """
    term1 = (
        viscosity_layer1 * (sinh(2.0 * theta1) - 2.0 * theta1) / viscosity_layer2
        / (cosh(2.0 * theta1) - 1.0 - 2.0 * theta1^2.0)
    )
    term2 = (
        (sinh(2.0 * theta2) - 2.0 * theta2)
        / (cosh(2.0 * theta2) - 1 - 2.0 * theta2^2.0)
    )
    d12_term = term1 + term2
    return d12_term
end

function calculate_j21_term(
    viscosity_layer1::Float64,
    viscosity_layer2::Float64,
    theta1::Float64,
    theta2::Float64
)::Float64
    """ Calculate j21 term for Rayleigh-Taylor instability benchmark.
    """
    term1 = (
        viscosity_layer1 * theta2 * (sinh(2.0 * theta1) + 2.0 * theta1)
        / viscosity_layer2 / (cosh(2.0 * theta1) - 1.0 - 2.0 * theta1^2.0)
    )
    term2 = (
        theta2 * (sinh(2.0 * theta2) + 2.0 * theta2)
        / (cosh(2.0 * theta2) - 1.0 - 2.0 * theta2^2.0)
    )
    j21_term = term1 + term2
    return j21_term
end

function calculate_j22_term(
    viscosity_layer1::Float64,
    viscosity_layer2::Float64,
    theta1::Float64,
    theta2::Float64
)::Float64
    """ Calculate j22 term for Rayleigh-Taylor instability benchmark.
    """
    term1 = (
        viscosity_layer1 * 2.0 * theta1^2.0 * theta2 / viscosity_layer2
        / (cosh(2.0 * theta1) - 1.0 - 2.0 * theta1^2.0)
    )
    term2 = 2.0 * theta2^3.0 / (cosh(2.0 * theta2) - 1.0 - 2.0 * theta2^2.0)
    j22_term = term1 - term2
    return j22_term
end

function print_benchmark_results(output_dict::OrderedDict{String, Any})::Nothing
    """ Print benchmark results.
    """
    println("")
    println(">> Rayleigh-Taylor Instability Benchmark Input and Output")
    for (key, value) in output_dict
        if occursin("cm/yr", key)
            println("   $(key): $(@sprintf("%.4E", value))")
        else
            println("   $(key): $(@sprintf("%.4f", value))")
        end
    end
    println("")
    return nothing
end

end # module 