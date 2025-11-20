module BoxConvectionManager

include("BoxConvection.jl")

import EarthBox.StaggeredGrid: Spacing
import EarthBox.ConversionFuncs: get_factor_cm_yr_to_m_s
import EarthBox.PlotToolsManager.Charts: plot_ncurves
import EarthBox.PlotToolsManager.Charts: make_plot_name
import EarthBox.EarthBoxDtypes: TestInfoDictType
import EarthBox.EarthBoxDtypes: ModelInputDictType
import EarthBox.EarthBoxDtypes: MaterialsDictType
import EarthBox.InputTools.Reader: make_parameters_dict
import ...Reader
import ...TestResults
import ...BenchmarksStruct: Benchmarks
import .BoxConvection

function compare_numerical_to_numerical(
    bench::Benchmarks
)::Tuple{Vector{Union{String, Float64}}, String}
    box_convection = BoxConvectionIsoviscous1a(bench)
    make_benchmark_plot(box_convection)
    return get_test_results(box_convection)
end

mutable struct BoxConvectionIsoviscous1a
    itime_step::Int
    main_paths::Dict{String, String}
    test_info_dict::Dict{String, Any}
    model_dict::ModelInputDictType
    materials_dict::Dict{Any, Any}
    materials_library_dict::MaterialsDictType
    tmyr::Float64
    relative_error_limit_percentage::Float64
    numerical_array::Vector{Float64}
    benchmark_array::Vector{Float64}
    benchmark_tmyr::Vector{Float64}
    tmyr_vec::Vector{Float64}
    rms_numerical::Vector{Float64}
    temp_avg_numerical::Vector{Float64}
    nusselt_numerical::Vector{Float64}
    rms_benchmark::Vector{Float64}
    nusselt_benchmark::Vector{Float64}
end

function BoxConvectionIsoviscous1a(bench::Benchmarks)
    relative_error_limit_percentage = 0.2
    numerical_array = zeros(Float64, 1)
    benchmark_array = zeros(Float64, 1)
    nsteps = bench.itime_step - 1
    tmyr_vec = zeros(Float64, nsteps)
    rms_numerical = zeros(Float64, nsteps)
    temp_avg_numerical = zeros(Float64, nsteps)
    nusselt_numerical = zeros(Float64, nsteps)
    rms_benchmark = zeros(Float64, 2)
    nusselt_benchmark = zeros(Float64, 2)

    time_steps = Vector{Int}()
    itime_step_start = 2
    for i in 1:nsteps
        itime_step = itime_step_start + i - 1
        push!(time_steps, itime_step)
    end

    rms_benchmark[1] = 42.865
    rms_benchmark[2] = 42.865

    nusselt_benchmark[1] = 4.884
    nusselt_benchmark[2] = 4.884

    icount = 1
    for itime_step in time_steps
        input_dict = get_inputs(itime_step, bench)
        tmyr_vec[icount] = input_dict["tmyr"]
        (
            rms_numerical[icount],
            temp_avg_numerical[icount],
            nusselt_numerical[icount]
        ) = BoxConvection.box_convection_quantities(input_dict)
        icount += 1
    end

    tmyr_final = tmyr_vec[nsteps]
    rms_numerical_final = rms_numerical[nsteps]

    benchmark_tmyr = [0.0, tmyr_final]

    numerical_array[1] = rms_numerical_final
    benchmark_array[1] = rms_benchmark[1]

    return BoxConvectionIsoviscous1a(
        bench.itime_step,
        bench.main_paths,
        bench.test_info_dict,
        bench.model_dict,
        bench.materials_dict,
        bench.materials_library_dict,
        tmyr_final,
        relative_error_limit_percentage,
        numerical_array,
        benchmark_array,
        benchmark_tmyr,
        tmyr_vec,
        rms_numerical,
        temp_avg_numerical,
        nusselt_numerical,
        rms_benchmark,
        nusselt_benchmark
    )
end

function get_inputs(
    itime_step::Int,
    bench::Benchmarks
)::Dict{String, Any}
    file_base_name = "vel_cmyr"
    (
        gridx_b_km,
        gridy_b_km,
        velocity_x_cm_yr,
        velocity_y_cm_yr,
        tmyr
    ) = Reader.read_vector_grid(
        bench.main_paths["post_proc_input_path"],
        file_base_name,
        itime_step
    )

    cm_yr_to_m_s = get_factor_cm_yr_to_m_s()
    vx1 = velocity_x_cm_yr .* cm_yr_to_m_s
    vy1 = velocity_y_cm_yr .* cm_yr_to_m_s

    gridx_b = gridx_b_km .* 1000.0
    gridy_b = gridy_b_km .* 1000.0
    xnum = length(gridx_b_km)
    ynum = length(gridy_b_km)
    xsize = gridx_b[xnum]
    ysize = gridy_b[ynum]

    xstp_b, ystp_b = Spacing.initialize_spacing_arrays(gridx_b, gridy_b)
    Spacing.update_basic_1Dgrid_spacing!(xnum, gridx_b, xstp_b)
    Spacing.update_basic_1Dgrid_spacing!(ynum, gridy_b, ystp_b)

    file_base_name = "TempC"
    (
        gridx_b_km,
        gridy_b_km,
        temperature_celsius,
        tmyr
    ) = Reader.read_scalar_grid(
        bench.main_paths["post_proc_input_path"],
        file_base_name,
        itime_step
    )

    tk1 = temperature_celsius .+ 273.0

    parameters_dict = make_parameters_dict(bench.model_dict)
    temperature_top_celsius = parameters_dict["temperature_top"][1]
    temperature_top = temperature_top_celsius + 273.0
    temperature_bottom_celsius = parameters_dict["temperature_bottom"][1]
    temperature_bottom = temperature_bottom_celsius + 273.0

    matid = 1
    material_name = bench.materials_dict[matid]["mat_name"][1]
    material_dict_layer1 = bench.materials_library_dict[material_name]
    density = material_dict_layer1["standard_density"]
    heat_capacity = material_dict_layer1["heat_capacity"]
    conductivity = material_dict_layer1["thermal_conductivity_ref"]

    return Dict{String, Any}(
        "tmyr" => tmyr,
        "gridx_b" => gridx_b,
        "gridy_b" => gridy_b,
        "xstp_b" => xstp_b,
        "ystp_b" => ystp_b,
        "xnum" => xnum,
        "ynum" => ynum,
        "xsize" => xsize,
        "ysize" => ysize,
        "vx1" => vx1,
        "vy1" => vy1,
        "tk1" => tk1,
        "temperature_top" => temperature_top,
        "temperature_bottom" => temperature_bottom,
        "density" => density,
        "heat_capacity" => heat_capacity,
        "conductivity" => conductivity
    )
end

function make_benchmark_plot(data::BoxConvectionIsoviscous1a)::Nothing
    data_xy = Dict{String, Any}(
        "x_arrays" => [
            data.benchmark_tmyr,
            data.tmyr_vec
        ],
        "y_arrays" => [
            data.rms_benchmark,
            data.rms_numerical
        ],
        "labels" => [
            "Benchmark",
            "Numerical"
        ],
        "colors" => [:blue, :red],
        "line_styles" => [:dash, :solid],
        "line_widths" => [2, 2],
        "line_colors" => [:blue, :red],
        "marker_sizes" => [10, 5],
        "marker_edge_colors" => [:black, :black],
        "marker_edge_widths" => [0.5, 0.5],
        "fill_styles" => [:none, :none]
    )

    axis_labels = ["Time (Myr)", "rms velocity"]
    plot_dimensions_xy = [0.0, 5000.0, 0.0, 70.0]

    boxtext_x = plot_dimensions_xy[2] * 0.5
    boxtext_y = plot_dimensions_xy[4] * 0.05
    boxtext = ""
    boxtext_info = [boxtext_x, boxtext_y, boxtext]

    title_base_text = "Isoviscous Convection Benchmark 1a"
    title = "$(title_base_text)"

    plot_name = data.main_paths["model_name"]
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
        "legend_location" => :bottomright,
        "legendfontsize" => 12,
        "figsize" => (10, 5),
        "guidefontsize" => 12,
        "titlefontsize" => 15,
        "tickfontsize" => 10,
        "xtick_size" => 500.0,
        "ytick_size" => 5.0
    )

    plot_ncurves(chart_input)
    return nothing
end

function get_test_results(
    data::BoxConvectionIsoviscous1a
)::Tuple{Vector{Union{String, Float64}}, String}
    result, result_msg = TestResults.get_test_results_numerical_vs_analytical(
        data.numerical_array,
        data.benchmark_array,
        data.relative_error_limit_percentage / 100.0
    )
    return result, result_msg
end

end # module 