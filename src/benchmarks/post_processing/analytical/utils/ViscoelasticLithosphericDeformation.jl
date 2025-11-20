module ViscoelasticLithosphericDeformation

using Printf
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelStructureManager: ModelStructure
import EarthBox.ModelStructureManager: calculate_structure!
import EarthBox.ModelStructureManager: get_basic_structure
import EarthBox.ModelStructureManager.AverageDensity: calculate_average_density_of_layers
import EarthBox.ModelStructureManager.LocalIsostasy: calc_local_isostatic_topography_for_2dgrid
import EarthBox.PlotToolsManager.Charts: plot_ncurves, make_plot_name
import EarthBox: MathTools
import EarthBox.EarthBoxDtypes
import EarthBox.InputTools.Reader: make_parameters_dict
import ...BenchmarkTools
import ...TestResults
import ...BenchmarksStruct: Benchmarks
import ...Reader

function execute_postprocessing_steps(
    bench::Benchmarks,
    plot_title::String,
    topo_plot_y_dimensions::Tuple{Float64, Float64}
)::Tuple{Vector{Union{String, Float64}}, String}
    input_dict = get_inputs(bench)

    sticky_ids = Int16[1, 2]
    sticky_water_ids = Int16[2]
    rock_ids = Int16[3, 4, 5, 6, 7, 8, 9, 10, 11]
    crust_ids = Int16[3, 4, 9, 10]
    lithosphere_ids = Int16[5, 6, 7, 11]
    asthenosphere_ids = Int16[8]

    model_structure = ModelStructure(
        input_dict["marknum"],
        input_dict["marker_matid"],
        input_dict["marker_x_m"],
        input_dict["marker_y_m"],
        input_dict["xnum"],
        input_dict["gridx_b"],
        input_dict["xstp_avg"],
        input_dict["ysize"],
        sticky_ids,
        sticky_water_ids,
        rock_ids,
        crust_ids,
        lithosphere_ids,
        asthenosphere_ids
    )

    calculate_structure!(model_structure)

    (
        ytopo_numerical,
        ymoho_numerical,
        _ytopo_smooth_numerical,
        ylith_base_numerical,
        _thick_water_numerical,
        thick_crust_numerical,
        _thick_mantle_lith_numerical,
        _thick_asthenosphere_numerical,
        _thick_sediments_numerical,
    ) = get_basic_structure(model_structure)

    (
        rho_sticky_x,
        rho_crust_x,
        rho_lithosphere_x,
        rho_asthenosphere_x
    ) = calculate_average_density_of_layers(
        input_dict["xnum"],
        input_dict["ynum"],
        input_dict["gridy_b"],
        ytopo_numerical,
        ymoho_numerical,
        ylith_base_numerical,
        input_dict["rho_kg_m3"]
    )

    # X-index reference used to define the reference column for the local
    # isostatic equation
    x_fraction_isostatic_reference = 0.05
    xindex_reference = floor(Int, input_dict["xnum"] * x_fraction_isostatic_reference) + 1
    if xindex_reference >= input_dict["xnum"]
        xindex_reference = input_dict["xnum"]
    end
    x_reference = input_dict["gridx_b"][xindex_reference]

    rho_sticky, rho_crust, rho_asthenosphere = get_density_model(bench)

    println("Density of sticky kg/m3 (model): ", rho_sticky)
    println("Density of crust kg/m3 (model): ", rho_crust)
    println("Density of asthenosphere kg/m3 (model): ", rho_asthenosphere)

    crustal_thickness_reference = get_initial_crustal_thickness(bench)
    println("Crustal thickness reference: ", crustal_thickness_reference)

    x_fraction_thinned = 0.5
    xindex_thinned = floor(Int, input_dict["xnum"] * x_fraction_thinned) + 1
    crustal_thickness_thinned = thick_crust_numerical[xindex_thinned]
    println("Crustal thickness thinned: ", crustal_thickness_thinned)

    isostatic_topography_analytical = calculate_local_isostatic_water_depth(
        rho_asthenosphere, rho_crust, rho_sticky,
        crustal_thickness_reference, crustal_thickness_thinned
    )

    gridy_b = input_dict["gridy_b"]
    ynum = input_dict["ynum"]
    y_bottom = gridy_b[ynum]
    ytopo_local = calc_local_isostatic_topography_for_2dgrid(
        xindex_reference,
        x_reference,
        y_bottom,
        input_dict["xnum"],
        ytopo_numerical,
        ymoho_numerical,
        ylith_base_numerical,
        rho_sticky_x,
        rho_crust_x,
        rho_lithosphere_x,
        rho_asthenosphere_x
    )

    ymin_plot_topo = topo_plot_y_dimensions[1]
    ymax_plot_topo = topo_plot_y_dimensions[2]
    topo_reference = ytopo_numerical[xindex_reference]
    x_reference = input_dict["gridx_b"][xindex_reference]
    make_benchmark_topo_plot(
        bench.main_paths,
        x_reference,
        topo_reference,
        bench.itime_step,
        input_dict["tmyr"],
        input_dict["gridx_b"],
        ytopo_numerical,
        ytopo_local,
        input_dict["ytopo_marker_chain"],
        ymin_plot_topo,
        ymax_plot_topo,
        isostatic_topography_analytical,
        plot_title
    )

    make_benchmark_crust_plot(
        bench.main_paths,
        bench.itime_step,
        input_dict["tmyr"],
        input_dict["gridx_b"],
        thick_crust_numerical
    )

    ytopo_target = fill(isostatic_topography_analytical, 1)
    xindex = floor(Int, length(input_dict["gridx_b"]) * 0.5) + 1
    ytopo_model = fill(ytopo_numerical[xindex] - topo_reference, 1)

    println("")
    println("ytopo_analytical (local isostatic): ", ytopo_target)
    println("ytopo_numerical: ", ytopo_model)

    relative_error_limit_percentage = 3.0
    result, result_msg = TestResults.get_test_results_numerical_vs_analytical(
        ytopo_model,
        ytopo_target,
        relative_error_limit_percentage / 100.0
    )
    println(result_msg)
    println(result)
    return result, result_msg
end

"""
Get material properties for all materials.
"""
function get_density_model(bench::Benchmarks)::Tuple{Float64, Float64, Float64}
    matid = 1
    rho_sticky = get_density(matid, bench)

    matid = 3
    rho_crust = get_density(matid, bench)

    matid = 8
    rho_asthenosphere = get_density(matid, bench)
    return rho_sticky, rho_crust, rho_asthenosphere
end

function get_density(
    matid::Int,
    bench::Benchmarks
)::Float64
    material_name = bench.materials_dict[matid]["mat_name"][1]
    material_dict_individual = bench.materials_library_dict[material_name]
    density = material_dict_individual["standard_density"]
    return density
end

function get_initial_crustal_thickness(
    bench::Benchmarks
)::Float64
    parameters_dict = make_parameters_dict(bench.model_dict)
    thick_upper_crust = parameters_dict["thick_upper_crust"][1]
    thick_lower_crust = parameters_dict["thick_lower_crust"][1]
    thick_crust_initial = thick_upper_crust + thick_lower_crust
    return thick_crust_initial
end

""" Calculate isostatic water depth.
"""
function calculate_local_isostatic_water_depth(
    rho_mantle::Float64,
    rho_crust::Float64,
    rho_water::Float64,
    thickness_crust_reference::Float64,
    thickness_crust::Float64
)::Float64
    water_depth_analytical = (
        (rho_mantle - rho_crust)
        / (rho_mantle - rho_water)
        * (thickness_crust_reference - thickness_crust)
    )
    return water_depth_analytical
end

function get_xindex_min_ytopo(ytopo_numerical::Vector{Float64})::Int
    min_ytopo = 1e38
    min_ytopo_index = -1
    for i in 1:length(ytopo_numerical)
        if ytopo_numerical[i] < min_ytopo
            min_ytopo_index = i
            min_ytopo = ytopo_numerical[i]
        end
    end
    return min_ytopo_index
end

function get_inputs(bench::Benchmarks)::Dict{String, Any}
    file_base_name = "particles"
    (
        marker_x_m,
        marker_y_m,
        marker_matid,
        tmyr
    ) = Reader.read_marker_material_ids(
        bench.main_paths["post_proc_input_path"], file_base_name, bench.itime_step)

    file_base_name = "rho_kg_m3"
    (
        gridx_b_km,
        gridy_b_km,
        rho_kg_m3,
        tmyr
    ) = Reader.read_scalar_grid(
        bench.main_paths["post_proc_input_path"], file_base_name, bench.itime_step)

    file_base_name = "topo"
    (
        xtopo_marker_chain_full,
        ytopo_marker_chain_full
    ) = Reader.read_topography(
        bench.main_paths["post_proc_input_path"], file_base_name, bench.itime_step)
    println("bench.itime_step: ", bench.itime_step)
    # Print min/max values of topography chain
    @printf("xtopo_marker_chain_full min: %.2f max: %.2f\n", 
            minimum(xtopo_marker_chain_full), maximum(xtopo_marker_chain_full))
    @printf("ytopo_marker_chain_full min: %.2f max: %.2f\n",
            minimum(ytopo_marker_chain_full), maximum(ytopo_marker_chain_full))

    file_base_name = "topo"
    (
        xtopo_marker_chain_ini_full,
        ytopo_marker_chain_ini_full
    ) = Reader.read_topography(
        bench.main_paths["post_proc_input_path"], file_base_name, 1)

    gridx_b = copy(gridx_b_km) .* 1000.0
    gridy_b = copy(gridy_b_km) .* 1000.0
    xnum = length(gridx_b_km)
    ynum = length(gridy_b_km)
    xsize = gridx_b[xnum]
    ysize = gridy_b[ynum]
    marknum = length(marker_x_m)
    xstp_avg = xsize / (xnum - 1)

    ytopo_marker_chain = zeros(xnum)
    MathTools.linear_interp_vals!(
        xtopo_marker_chain_full,
        ytopo_marker_chain_full,
        gridx_b,
        ytopo_marker_chain
    )

    ytopo_marker_chain_ini = zeros(xnum)
    MathTools.linear_interp_vals!(
        xtopo_marker_chain_ini_full,
        ytopo_marker_chain_ini_full,
        gridx_b,
        ytopo_marker_chain_ini
    )

    input_dict = Dict{String, Any}(
        "tmyr" => tmyr,
        "gridx_b" => gridx_b,
        "gridy_b" => gridy_b,
        "rho_kg_m3" => rho_kg_m3,
        "xnum" => xnum,
        "ynum" => ynum,
        "xsize" => xsize,
        "ysize" => ysize,
        "xstp_avg" => xstp_avg,
        "marknum" => marknum,
        "marker_x_m" => marker_x_m,
        "marker_y_m" => marker_y_m,
        "marker_matid" => marker_matid,
        "ytopo_marker_chain" => ytopo_marker_chain,
        "ytopo_marker_chain_ini" => ytopo_marker_chain_ini
    )

    return input_dict
end

""" Calculate the average difference in local and numerical topography
"""
function calculate_average_central_topo_difference(
    water_depth_analytical::Float64,
    ytopo_numerical::Vector{Float64},
    topo_reference::Float64,
    gridx_b::Vector{Float64},
    width::Float64,
    x_fraction::Float64
)::Float64
    xindex = floor(Int, length(gridx_b) * x_fraction)
    xmid = gridx_b[xindex]
    xmin = xmid - width
    xmax = xmid + width
    xnum = length(ytopo_numerical)

    icount = 0
    sumit = 0.0
    for i in 1:xnum
        x = gridx_b[i]
        ytopo_numerical_val = ytopo_numerical[i] - topo_reference
        diff = abs(ytopo_numerical_val - water_depth_analytical)
        if xmin < x < xmax
            icount += 1
            sumit = sumit + diff
        end
    end
    diff_avg = sumit / icount
    return diff_avg
end

""" Custom function for making benchmark plots.
"""
function make_benchmark_topo_plot(
    main_paths::Dict{String, String},
    x_reference::Float64,
    topo_reference::Float64,
    itime_step::Int64,
    tmyr::Float64,
    gridx_b::Vector{Float64},
    ytopo_numerical::Vector{Float64},
    ytopo_local::Vector{Float64},
    ytopo_marker_chain::Vector{Float64},
    ymin_plot::Float64,
    ymax_plot::Float64,
    water_depth_analytical::Float64,
    plot_title::String
)::Nothing
    data_xy = Dict{String, Any}(
        "x_arrays" => [
            gridx_b ./ 1000.0,
            gridx_b ./ 1000.0,
            gridx_b ./ 1000.0,
            [x_reference / 1000.0],
            [0, gridx_b[end] / 1000.0]
        ],
        "y_arrays" => [
            ytopo_marker_chain ./ 1000.0 .- topo_reference / 1000.0,
            ytopo_numerical ./ 1000.0 .- topo_reference / 1000.0,
            ytopo_local ./ 1000.0 .- topo_reference / 1000.0,
            [0.0],
            [water_depth_analytical / 1000.0, water_depth_analytical / 1000.0]
        ],
        "labels" => [
            "Topography From Marker Chain Interpolated to X-Grid",
            "Topography From Markers",
            "Local Isostatic Topography Calculated From Markers",
            "Reference Marker Column Location for Isostatic Calculation",
            "Analytical Local Isostatic Solution"
        ],
        "colors" => [:black, :green, :red, :blue, :black],
        "line_styles" => [:solid, :solid, :solid, :solid, :dash],
        "line_widths" => [2, 1, 1, 1, 1],
        "line_colors" => [:black, :green, :red, :blue, :black],
        "marker_sizes" => [5, 5, 5, 5, 5],
        "marker_edge_colors" => [:black, :black, :black, :black, :black],
        "marker_edge_widths" => [0.5, 0.5, 0.5, 0.5, 0.5],
        "fill_styles" => [:none, :none, :none, :circle, :none]
    )

    xmax = gridx_b[end] / 1000.0
    axis_labels = ["x (km)", "Topography Relative to Reference (km)"]
    plot_dimensions_xy = [0, xmax, ymin_plot, ymax_plot]

    boxtext_x = plot_dimensions_xy[2] * 0.65
    boxtext_y = plot_dimensions_xy[4] * 0.98
    boxtext = ""
    boxtext_info = [boxtext_x, boxtext_y, boxtext]

    title_base_text = plot_title
    title = "$(title_base_text): Numerical Time $(@sprintf("%.2f", tmyr)) Myr"

    plot_name = make_plot_name(
        main_paths["model_name"] * "_topo",
        itime_step
    )
    plot_file_path = joinpath(
        main_paths["post_proc_output_path"], plot_name)

    chart_input = Dict{String, Any}(
        "plot_file_path" => plot_file_path,
        "title" => title,
        "axis_labels" => axis_labels,
        "plot_dimensions_xy" => plot_dimensions_xy,
        "data_xy" => data_xy,
        "boxtext_info" => boxtext_info,
        "iuse_inversion" => 1,
        "aspect_ratio" => :auto,
        "figure_dpi" => 150,
        "legend_location" => :bottomright,
        "legendfontsize" => 12,
        "figsize" => (10, 5),
        "guidefontsize" => 12,
        "titlefontsize" => 15,
        "tickfontsize" => 10,
        "annotation_fontsize" => 10,
        "xtick_size" => 50.0,
        "ytick_size" => 2.0
    )

    plot_ncurves(chart_input)
    return nothing
end

""" Custom function for making benchmark plots.
"""
function make_benchmark_crust_plot(
    main_paths::Dict{String, String},
    itime_step::Int,
    tmyr::Float64,
    gridx_b::Vector{Float64},
    thick_crust_numerical::Vector{Float64}
)::Nothing
    data_xy = Dict{String, Any}(
        "x_arrays" => [
            gridx_b ./ 1000.0,
        ],
        "y_arrays" => [
            thick_crust_numerical ./ 1000.0
        ],
        "labels" => [
            "Crustal Thickness From Markers",
        ],
        "colors" => [:green],
        "line_styles" => [:solid],
        "line_widths" => [1],
        "line_colors" => [:green],
        "marker_sizes" => [5],
        "marker_edge_colors" => [:black],
        "marker_edge_widths" => [0.5],
        "fill_styles" => [:none]
    )

    axis_labels = [
        "x (km)",
        "Crustal Thickness (km)"
    ]
    xmax = gridx_b[end] / 1000.0
    plot_dimensions_xy = [0, xmax, 5, 35.0]

    boxtext_x = plot_dimensions_xy[2] * 0.65
    boxtext_y = plot_dimensions_xy[4] * 0.98
    boxtext = ""
    boxtext_info = [boxtext_x, boxtext_y, boxtext]

    title_base_text = "Viscoelastic Extension Benchmark"
    title = "$(title_base_text): Numerical Time $(@sprintf("%.2f", tmyr)) Myr"

    plot_name = make_plot_name(
        main_paths["model_name"] * "_xth",
        itime_step
    )
    plot_file_path = joinpath(
        main_paths["post_proc_output_path"], plot_name)

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
        "figsize" => (10, 5),
        "guidefontsize" => 12,
        "titlefontsize" => 15,
        "tickfontsize" => 10,
        "annotation_fontsize" => 10,
        "xtick_size" => 50.0,
        "ytick_size" => 2.0
    )

    plot_ncurves(chart_input)
    return nothing
end

end # module 