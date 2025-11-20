module SeafloorSpreadingManager

import Printf: @sprintf
import EarthBox.ModelStructureManager: ModelStructure
import EarthBox.ModelStructureManager: calculate_structure!
import EarthBox.ModelStructureManager: get_basic_structure
import EarthBox.ModelStructureManager.AverageDensity: calculate_average_density_of_layers
import EarthBox.ModelStructureManager.LocalIsostasy: calc_local_isostatic_topography_for_2dgrid
import EarthBox.PlotToolsManager.Charts: plot_ncurves
import EarthBox.PlotToolsManager.Charts: make_plot_name
import EarthBox.EarthBoxDtypes: TestInfoDictType
import EarthBox.EarthBoxDtypes: ModelInputDictType
import EarthBox.EarthBoxDtypes: MaterialsDictType
import EarthBox.InputTools.Reader: make_parameters_dict
import EarthBox.MathTools: linear_interp_vals!
import ...Reader
import ...TestResults
import ...BenchmarksStruct: Benchmarks

function compare_numerical_to_empirical(
    bench::Benchmarks
)::Tuple{Vector{Union{String, Float64}}, String}
    input_dict = get_inputs(bench.itime_step, bench)

    model_structure = ModelStructure(
        input_dict["marknum"],
        input_dict["marker_matid"],
        input_dict["marker_x_m"],
        input_dict["marker_y_m"],
        input_dict["xnum"],
        input_dict["gridx_b"],
        input_dict["xstp_avg"],
        input_dict["ysize"],
        Int16[1, 2], # matids_sticky
        Int16[2], # matids_sticky_water
        Int16[3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], # matids_rock
        Int16[8, 9, 10, 11, 12, 13, 14, 15], # matids_crust
        Int16[4], # matids_lithosphere
        Int16[5, 6, 7], # matids_asthenosphere
        Int16[-1] # matids_sediments
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
        _thick_sediments_numerical
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

    # Get the x-index of minimum topographic depth to approximate the location
    # of the ridge axis.
    xindex_min_ytopo = get_xindex_min_ytopo(ytopo_numerical)
    x_min_ytopo = input_dict["gridx_b"][xindex_min_ytopo]
    println(
        "xindex_max_topo: ", xindex_min_ytopo,
        " x_min_ytopo (km): ", x_min_ytopo/1000.0
    )

    # X-index reference used to define the reference column for the local
    # isostatic equation
    xindex_reference = floor(Int, xindex_min_ytopo/4)
    if xindex_reference >= input_dict["xnum"]
        xindex_reference = input_dict["xnum"] - 1
    end
    x_reference = input_dict["gridx_b"][xindex_reference]

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

    make_benchmark_topo_plot(
        bench.main_paths,
        xindex_reference,
        bench.itime_step,
        input_dict["tmyr"],
        input_dict["gridx_b"],
        ytopo_numerical,
        ytopo_local,
        input_dict["ytopo_marker_chain"]
    )

    make_benchmark_crust_plot(
        bench.main_paths,
        bench.itime_step,
        input_dict["tmyr"],
        input_dict["gridx_b"],
        thick_crust_numerical
    )

    thick_oc_numerical_avg_value = calculate_average_thickness_excluding_ridge(
        thick_crust_numerical,
        input_dict["gridx_b"],
        x_min_ytopo
    )

    thick_oc_numerical_avg = [thick_oc_numerical_avg_value]
    thick_oc_avg_target = [7000.0]

    println("")
    println(
        "Average off-axis thickness of oceanic crust (Numerical): ",
        thick_oc_numerical_avg[1]
    )
    println(
        "Average thickness of oceanic crust (Target): ",
        thick_oc_avg_target[1]
    )

    relative_error_limit_percentage = 6.0
    result, result_msg = TestResults.get_test_results_numerical_vs_analytical(
        thick_oc_numerical_avg,
        thick_oc_avg_target,
        relative_error_limit_percentage/100.0
    )
    return result, result_msg
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

function calculate_average_thickness_excluding_ridge(
    thick_crust_numerical::Vector{Float64},
    gridx_b::Vector{Float64},
    x_min_ytopo::Float64
)::Float64
    ridge_width = 25000.0
    xmin_ridge = x_min_ytopo - ridge_width
    xmax_ridge = x_min_ytopo + ridge_width
    xnum = length(thick_crust_numerical)

    icount = 0
    sumit = 0.0
    for i in 1:xnum
        x = gridx_b[i]
        thick = thick_crust_numerical[i]
        if x < xmin_ridge || x > xmax_ridge
            icount += 1
            sumit += thick
        end
    end
    thick_crust_numerical_avg = sumit/icount
    return thick_crust_numerical_avg
end

function get_inputs(
    itime_step::Int,
    bench::Benchmarks
)::Dict{String, Any}
    file_base_name = "particles"
    (
        marker_x_m,
        marker_y_m,
        marker_matid,
        tmyr
    ) = Reader.read_marker_material_ids(
        bench.main_paths["post_proc_input_path"],
        file_base_name,
        itime_step
    )

    file_base_name = "rho_kg_m3"
    (
        gridx_b_km,
        gridy_b_km,
        rho_kg_m3,
        tmyr
    ) = Reader.read_scalar_grid(
        bench.main_paths["post_proc_input_path"],
        file_base_name,
        itime_step
    )

    file_base_name = "topo"
    (
        xtopo_marker_chain_full,
        ytopo_marker_chain_full
    ) = Reader.read_topography(
        bench.main_paths["post_proc_input_path"],
        file_base_name,
        itime_step
    )

    file_base_name = "topo"
    (
        xtopo_marker_chain_ini_full,
        ytopo_marker_chain_ini_full
    ) = Reader.read_topography(
        bench.main_paths["post_proc_input_path"],
        file_base_name,
        1
    )

    gridx_b = gridx_b_km .* 1000.0
    gridy_b = gridy_b_km .* 1000.0
    xnum = length(gridx_b_km)
    ynum = length(gridy_b_km)
    xsize = gridx_b[xnum]
    ysize = gridy_b[ynum]
    marknum = length(marker_x_m)
    xstp_avg = xsize/(xnum - 1)

    ytopo_marker_chain = zeros(Float64, xnum)
    linear_interp_vals!(
        xtopo_marker_chain_full,
        ytopo_marker_chain_full,
        gridx_b,
        ytopo_marker_chain
    )

    ytopo_marker_chain_ini = zeros(Float64, xnum)
    linear_interp_vals!(
        xtopo_marker_chain_ini_full,
        ytopo_marker_chain_ini_full,
        gridx_b,
        ytopo_marker_chain_ini
    )

    return Dict{String, Any}(
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
end

function make_benchmark_topo_plot(
    main_paths::Dict{String, String},
    xindex_reference::Int,
    itime_step::Int,
    tmyr::Float64,
    gridx_b::Vector{Float64},
    ytopo_numerical::Vector{Float64},
    ytopo_local::Vector{Float64},
    ytopo_marker_chain::Vector{Float64}
)::Nothing
    topo_ref = ytopo_numerical[xindex_reference]/1000.0
    x_reference = gridx_b[xindex_reference]/1000.0

    data_xy = Dict{String, Any}(
        "x_arrays" => [
            gridx_b ./ 1000.0,
            gridx_b ./ 1000.0,
            gridx_b ./ 1000.0,
            [x_reference]
        ],
        "y_arrays" => [
            ytopo_marker_chain ./ 1000.0 .- topo_ref,
            ytopo_numerical ./ 1000.0 .- topo_ref,
            ytopo_local ./ 1000.0 .- topo_ref,
            [0.0]
        ],
        "labels" => [
            "Topography Marker Chain Interpolated to X-Grid",
            "Topography From Markers",
            "Local Isostatic Topography From Markers",
            "Reference Marker Column Location"
        ],
        "colors" => [:black, :green, :red, :blue],
        "line_styles" => [:solid, :solid, :solid, :dot],
        "line_widths" => [2, 2, 2, 2],
        "line_colors" => [:black, :green, :red, :blue],
        "marker_sizes" => [5, 5, 5, 10],
        "marker_edge_colors" => [:black, :black, :black, :black],
        "marker_edge_widths" => [0.5, 0.5, 0.5, 0.5],
        "fill_styles" => [:none, :none, :none, :none]
    )

    axis_labels = ["x (km)", "Topography Relative to Reference (km)"]
    plot_dimensions_xy = [0.0, 500.0, -5.0, 5.0]

    boxtext_x = plot_dimensions_xy[2] * 0.65
    boxtext_y = plot_dimensions_xy[4] * 0.98
    boxtext = ""
    boxtext_info = [boxtext_x, boxtext_y, boxtext]

    title_base_text = "Seafloor Spreading Benchmark"
    title = "$(title_base_text): Numerical Time $(@sprintf("%.2f", tmyr)) Myr"

    plot_name = make_plot_name(
        main_paths["model_name"] * "_topo",
        itime_step
    )
    plot_file_path = joinpath(
        main_paths["post_proc_output_path"],
        plot_name
    )

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
        "xtick_size" => 50.0,
        "ytick_size" => 1.0
    )

    plot_ncurves(chart_input)
    return nothing
end

function make_benchmark_crust_plot(
    main_paths::Dict{String, String},
    itime_step::Int,
    tmyr::Float64,
    gridx_b::Vector{Float64},
    thick_crust_numerical::Vector{Float64}
)::Nothing
    data_xy = Dict{String, Any}(
        "x_arrays" => [
            gridx_b ./ 1000.0
        ],
        "y_arrays" => [
            thick_crust_numerical ./ 1000.0
        ],
        "labels" => [
            "Solidified and Partially Molten Gabbro + Gabbroic Magma"
        ],
        "colors" => [:green],
        "line_styles" => [:solid],
        "line_widths" => [2],
        "line_colors" => [:green],
        "marker_sizes" => [5],
        "marker_edge_colors" => [:black],
        "marker_edge_widths" => [0.5],
        "fill_styles" => [:none]
    )

    axis_labels = ["x (km)", "Crustal Thickness (km)"]
    plot_dimensions_xy = [0.0, 500.0, 1.0, 15.0]

    boxtext_x = plot_dimensions_xy[2] * 0.65
    boxtext_y = plot_dimensions_xy[4] * 0.98
    boxtext = ""
    boxtext_info = [boxtext_x, boxtext_y, boxtext]

    title_base_text = "Seafloor Spreading Benchmark"
    title = "$(title_base_text): Numerical Time $(@sprintf("%.2f", tmyr)) Myr"

    plot_name = make_plot_name(
        main_paths["model_name"] * "_xth",
        itime_step
    )
    plot_file_path = joinpath(
        main_paths["post_proc_output_path"],
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
        "legend_location" => :topleft,
        "legendfontsize" => 12,
        "figsize" => (10, 5),
        "guidefontsize" => 12,
        "titlefontsize" => 15,
        "tickfontsize" => 10,
        "xtick_size" => 50.0,
        "ytick_size" => 0.5
    )

    plot_ncurves(chart_input)
    return nothing
end

end # module 