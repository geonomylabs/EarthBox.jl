module FlexureTriangularHole

include("LithosFlexure.jl")

using Printf
import EarthBox.ModelStructureManager: ModelStructure
import EarthBox.ModelStructureManager: calculate_structure!
import EarthBox.ModelStructureManager: get_basic_structure
import EarthBox.PlotToolsManager.Charts: plot_ncurves, make_plot_name
import EarthBox: MathTools
import EarthBox: EarthBoxDtypes
import EarthBox.InputTools.Reader: make_parameters_dict
import .LithosFlexure: calc_flexural_deflection_from_block_loads
import ...BenchmarkTools
import ...TestResults
import ...BenchmarksStruct: Benchmarks
import ...Reader

function compare_numerical_to_analytical(
    bench::Benchmarks
)::Tuple{Vector{Union{String, Float64}}, String}
    flexure_triangular_hole = FlexureTriangularHoleStruct(bench)
    make_benchmark_plot(flexure_triangular_hole)
    return get_test_results(flexure_triangular_hole)
end

mutable struct FlexureTriangularHoleStruct
    itime_step::Int
    main_paths::Dict{String, String}
    test_info_dict::Dict{String, Any}
    model_dict::EarthBoxDtypes.ModelInputDictType
    materials_dict::Dict{Any, Any}
    materials_library_dict::EarthBoxDtypes.MaterialsDictType
    tmyr::Float64
    relative_error_limit_percentage::Float64
    input_dict::Dict{String, Any}
    ytopo_numerical::Vector{Float64}
    ymoho_numerical::Vector{Float64}
    ytopo_numerical_smooth::Vector{Float64}
    ylith_base_numerical::Vector{Float64}
    thick_crust_numerical::Vector{Float64}
    ytopo_regional_analytical::Vector{Float64}
    delta::Float64
end

function FlexureTriangularHoleStruct(bench::Benchmarks)
    input_dict = get_inputs(bench)

    structure_from_markers = get_structure_from_markers(input_dict)
    ytopo_numerical = structure_from_markers["ytopo"]
    ymoho_numerical = structure_from_markers["ymoho"]
    ytopo_numerical_smooth = structure_from_markers["ytopo_smooth"]
    ylith_base_numerical = structure_from_markers["ylith_base"]
    thick_crust_numerical = structure_from_markers["thick_crust"]

    add_structure_inputs!(input_dict, ytopo_numerical, ymoho_numerical, 
        ytopo_numerical_smooth, ylith_base_numerical, thick_crust_numerical)
    print_inputs(input_dict)

    (
        ytopo_regional_analytical, _xcoors_regional_analytical
    ) = analytical_flexure(input_dict)

    delta = maximum(abs.(
        input_dict["ytopo_marker_chain"] - input_dict["ytopo_numerical_smooth"]
    ))
    println("delta: ", delta)

    # Note that the python version had a max relative error of 8.0%
    return FlexureTriangularHoleStruct(
        bench.itime_step, bench.main_paths, bench.test_info_dict, bench.model_dict,
        bench.materials_dict, bench.materials_library_dict, input_dict["tmyr"],
        15.0, input_dict, ytopo_numerical,
        ymoho_numerical, ytopo_numerical_smooth, ylith_base_numerical,
        thick_crust_numerical, ytopo_regional_analytical, delta
    )
end

function get_inputs(bench::Benchmarks)::Dict{String, Any}
    itime_step = bench.itime_step
    main_paths = bench.main_paths
    model_dict = bench.model_dict
    materials_dict = bench.materials_dict
    materials_library_dict = bench.materials_library_dict

    poissons_ratio = 0.5

    t_crust_ini, grav_acceleration = get_model_info(model_dict)
    thick_elastic_km = t_crust_ini/1000.0

    (
        rho_load, rho_infill, rho_disp, rho_mantle, youngs_modulus
    ) = get_all_properties(materials_dict, materials_library_dict)

    file_base_name = "particles"
    (
        marker_x_m, marker_y_m, marker_matid, tmyr
    ) = Reader.read_marker_material_ids(
        main_paths["post_proc_input_path"], file_base_name, itime_step)

    file_base_name = "rho_kg_m3"
    (
        gridx_b_km, gridy_b_km, rho_kg_m3, tmyr
    ) = Reader.read_scalar_grid(
            main_paths["post_proc_input_path"], file_base_name, itime_step)

    gridx_b = copy(gridx_b_km) .* 1000.0
    gridy_b = copy(gridy_b_km) .* 1000.0
    xnum = length(gridx_b_km)
    ynum = length(gridy_b_km)
    xsize = gridx_b[xnum]
    ysize = gridy_b[ynum]
    marknum = length(marker_x_m)
    xstp_avg = xsize / (xnum - 1)

    (
        ytopo_marker_chain, ytopo_marker_chain_ini
    ) = get_topo_data_interpolated_to_grid_coordinates(
            xnum, gridx_b, main_paths, itime_step)

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
        "ytopo_marker_chain_ini" => ytopo_marker_chain_ini,
        "rho_mantle" => rho_mantle,
        "rho_load" => rho_load,
        "rho_disp" => rho_disp,
        "rho_infill" => rho_infill,
        "youngs_modulus" => youngs_modulus,
        "t_crust_ini" => t_crust_ini,
        "thick_elastic_km" => thick_elastic_km,
        "grav_acceleration" => grav_acceleration,
        "poissons_ratio" => poissons_ratio
    )
    return input_dict
end

function get_topo_data_interpolated_to_grid_coordinates(
    xnum::Int,
    gridx_b::Vector{Float64},
    main_paths::Dict{String, String},
    itime_step::Int
)::Tuple{Vector{Float64}, Vector{Float64}}
    file_base_name = "topo"
    (
        xtopo_marker_chain_full, ytopo_marker_chain_full
    ) = Reader.read_topography(
            main_paths["post_proc_input_path"], file_base_name, itime_step)

    file_base_name = "topo"
    (
        xtopo_marker_chain_ini_full, ytopo_marker_chain_ini_full
    ) = Reader.read_topography(
            main_paths["post_proc_input_path"], file_base_name, 1)

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
    return ytopo_marker_chain, ytopo_marker_chain_ini
end

function print_inputs(input_dict::Dict{String, Any})::Nothing
    println("grav_acceleration: ", input_dict["grav_acceleration"])
    println("t_crust_ini: ", input_dict["t_crust_ini"])
    println("thick_elastic_km: ", input_dict["thick_elastic_km"])
    println("rho_load: ", input_dict["rho_load"])
    println("rho_infill: ", input_dict["rho_infill"])
    println("rho_disp: ", input_dict["rho_disp"])
    println("rho_mantle: ", input_dict["rho_mantle"])
    println("youngs_modulus: ", input_dict["youngs_modulus"])
    println("poissons_ratio: ", input_dict["poissons_ratio"])
    ytopo_numerical = input_dict["ytopo_numerical"]
    ytopo_marker_chain_ini = input_dict["ytopo_marker_chain_ini"]
    thick_crust = input_dict["thick_crust_numerical"]
    println("ytopo_numerical min: ", minimum(ytopo_numerical))
    println("ytopo_numerical max: ", maximum(ytopo_numerical))
    println("ytopo_marker_chain_ini min: ", minimum(ytopo_marker_chain_ini))
    println("ytopo_marker_chain_ini max: ", maximum(ytopo_marker_chain_ini))
    println("thick_crust min: ", minimum(thick_crust))
    println("thick_crust max: ", maximum(thick_crust))
    println("ytopo_numerical[1]: ", ytopo_numerical[1])
    println("ytopo_marker_chain_ini[1]: ", ytopo_marker_chain_ini[1])
    return nothing
end

function get_all_properties(
    materials_dict::Union{Dict{Any, Any}, Nothing},
    materials_library_dict::EarthBoxDtypes.MaterialsDictType
)::Tuple{Float64, Float64, Float64, Float64, Float64}
    matid = 1
    rho_load, _shear_modulus = get_properties(matid, materials_dict, 
        materials_library_dict)
    rho_infill = rho_load
    matid = 4
    rho_disp, youngs_modulus = get_properties(matid, materials_dict, 
        materials_library_dict)
    matid = 6
    rho_mantle, _shear_modulus = get_properties(matid, materials_dict, 
        materials_library_dict)
    return rho_load, rho_infill, rho_disp, rho_mantle, youngs_modulus
end

function get_properties(
    matid::Int,
    materials_dict::Union{Dict{Any, Any}, Nothing},
    materials_library_dict::EarthBoxDtypes.MaterialsDictType
)::Tuple{Float64, Float64}
    material_name = materials_dict[matid]["mat_name"][1]
    material_dict_individual = materials_library_dict[material_name]
    density = material_dict_individual["standard_density"]
    shear_modulus = material_dict_individual["shear_modulus"]
    return density, shear_modulus
end

function get_model_info(
    model_dict::Union{EarthBoxDtypes.ModelInputDictType, Nothing}
)::Tuple{Float64, Float64}
    parameters_dict = make_parameters_dict(model_dict)
    thick_upper_crust = parameters_dict["thick_upper_crust"][1]
    thick_lower_crust = parameters_dict["thick_lower_crust"][1]
    thick_crust_initial = thick_upper_crust + thick_lower_crust
    grav_acceleration = parameters_dict["gravity_y"][1]
    return thick_crust_initial, grav_acceleration
end

function get_structure_from_markers(
    input_dict::Dict{String, Any}
)::Dict{String, Vector{Float64}}
    model_structure = ModelStructure(
        input_dict["marknum"],
        input_dict["marker_matid"],
        input_dict["marker_x_m"],
        input_dict["marker_y_m"],
        input_dict["xnum"],
        input_dict["gridx_b"],
        input_dict["xstp_avg"],
        input_dict["ysize"],
        Int16[1, 2],
        Int16[2],
        Int16[3, 4, 5, 6, 7, 8, 9],
        Int16[4, 5],
        Int16[6, 7, 8],
        Int16[9]
    )

    calculate_structure!(model_structure)

    (
        ytopo,
        ymoho,
        ytopo_smooth,
        ylith_base,
        _thick_water,
        thick_crust,
        _thick_mantle_lith,
        _thick_asthenosphere,
        _thick_seds
    ) = get_basic_structure(model_structure)

    structure_from_markers = Dict{String, Vector{Float64}}(
        "ytopo" => ytopo,
        "ymoho" => ymoho,
        "ytopo_smooth" => ytopo_smooth,
        "ylith_base" => ylith_base,
        "thick_crust" => thick_crust
    )

    return structure_from_markers
end


function add_structure_inputs!(
    input_dict::Dict{String, Any},
    ytopo_numerical::Vector{Float64},
    ymoho_numerical::Vector{Float64},
    ytopo_numerical_smooth::Vector{Float64},
    ylith_base_numerical::Vector{Float64},
    thick_crust_numerical::Vector{Float64}
)::Nothing
    input_dict["ytopo_numerical"] = ytopo_numerical
    input_dict["ymoho_numerical"] = ymoho_numerical
    input_dict["ytopo_numerical_smooth"] = ytopo_numerical_smooth
    input_dict["ylith_base_numerical"] = ylith_base_numerical
    input_dict["thick_crust_numerical"] = thick_crust_numerical
    return nothing
end

function make_benchmark_plot(data::FlexureTriangularHoleStruct)::Nothing
    tmyr = data.input_dict["tmyr"]
    gridx_b = data.input_dict["gridx_b"]
    ytopo_marker_chain = data.input_dict["ytopo_marker_chain"]

    data_xy = Dict{String, Any}(
        "x_arrays" => [
            gridx_b ./ 1000.0,
            gridx_b ./ 1000.0,
            gridx_b ./ 1000.0,
            gridx_b ./ 1000.0
        ],
        "y_arrays" => [
            data.ytopo_numerical ./ 1000.0,
            ytopo_marker_chain ./ 1000.0,
            data.ytopo_regional_analytical ./ 1000.0,
            data.ytopo_numerical_smooth ./ 1000.0
        ],
        "labels" => [
            "Numerical_Markers",
            "Numerical_Marker_Chain",
            "Analytical_Regional",
            "Numerical_Markers_LPF"
        ],
        "colors" => [:green, :red, :blue, :pink],
        "line_styles" => [:solid, :solid, :dash, :solid],
        "line_widths" => [2, 2, 2, 2],
        "line_colors" => [:green, :red, :blue, :pink],
        "marker_sizes" => [5, 5, 5, 5],
        "marker_edge_colors" => [:black, :black, :black, :black],
        "marker_edge_widths" => [0.5, 0.5, 0.5, 0.5],
        "fill_styles" => [:none, :none, :none, :none]
    )

    axis_labels = ["x (km)", "Y (km)"]
    plot_dimensions_xy = [0, 600.0, 5.0, 30.0]

    boxtext_x = plot_dimensions_xy[2] * 0.65
    boxtext_y = plot_dimensions_xy[4] * 0.98
    boxtext = ""
    boxtext_info = [boxtext_x, boxtext_y, boxtext]

    title_base_text = "Flexure Crustal Hole Load"
    title = "$(title_base_text): Numerical Time $(@sprintf("%.2f", tmyr)) Myr"

    plot_name = make_plot_name(
        data.main_paths["model_name"] * "_topo",
        data.itime_step
    )
    plot_file_path = joinpath(
        data.main_paths["post_proc_output_path"], plot_name)

    chart_input = Dict{String, Any}(
        "plot_file_path" => plot_file_path,
        "title" => title,
        "axis_labels" => axis_labels,
        "plot_dimensions_xy" => plot_dimensions_xy,
        "data_xy" => data_xy,
        "boxtext_info" => boxtext_info,
        "iuse_inversion" => 1,
        "aspect_ratio" => 10,
        "figure_dpi" => 150,
        "legend_location" => :bottomleft,
        "legendfontsize" => 12,
        "figsize" => (10, 5),
        "guidefontsize" => 12,
        "titlefontsize" => 15,
        "tickfontsize" => 10,
        "annotation_fontsize" => 10,
        "xtick_size" => 50.0,
        "ytick_size" => 5.0
    )

    plot_ncurves(chart_input)
    return nothing
end

function get_test_results(
    data::FlexureTriangularHoleStruct
)::Tuple{Vector{Union{String, Float64}}, String}
    result, result_msg = TestResults.get_test_results_numerical_vs_analytical(
        data.ytopo_numerical,
        data.ytopo_regional_analytical,
        data.relative_error_limit_percentage/100.0
    )
    return result, result_msg
end

function analytical_flexure(
    input_dict::Dict{String, Any}
)::Tuple{Vector{Float64}, Vector{Float64}}
    t_crust_ini = input_dict["t_crust_ini"]
    thick_elastic_km = input_dict["thick_elastic_km"]
    youngs_modulus = input_dict["youngs_modulus"]
    poissons_ratio = input_dict["poissons_ratio"]
    grav_acceleration = input_dict["grav_acceleration"]
    rho_mantle = input_dict["rho_mantle"]
    rho_load = input_dict["rho_load"]
    rho_disp = input_dict["rho_disp"]
    rho_infill = input_dict["rho_infill"]

    xnum = input_dict["xnum"]
    gridx_b = input_dict["gridx_b"]
    ytopo = input_dict["ytopo_numerical"]
    ytopo_marker_chain_ini = input_dict["ytopo_marker_chain_ini"]
    thick_crust = input_dict["thick_crust_numerical"]

    xcoor_blocks = zeros(xnum)
    width_blocks = zeros(xnum)
    height_blocks = zeros(xnum)
    for j in 1:xnum
        xcoor_blocks[j] = gridx_b[j]
        height_blocks[j] = t_crust_ini - thick_crust[j]
        if j == 1
            width_blocks[j] = (gridx_b[j+1] - gridx_b[j])/2.0
        elseif 1 < j < xnum
            width_blocks[j] = (
                (gridx_b[j] - gridx_b[j-1])/2.0
                + (gridx_b[j+1] - gridx_b[j])/2.0
            )
        elseif j == xnum
            width_blocks[j] = gridx_b[j] - gridx_b[j-1]
        end
    end

    deflections = calc_flexural_deflection_from_block_loads(
        xcoor_blocks,
        width_blocks,
        height_blocks,
        thick_elastic_km,
        youngs_modulus,
        poissons_ratio,
        grav_acceleration,
        rho_load,
        rho_disp,
        rho_mantle,
        rho_infill
    )

    println("deflections[1]: ", deflections[1])
    println("ytopo_marker_chain_ini[1]: ", ytopo_marker_chain_ini[1])
    println("ytopo[1] - ytopo_marker_chain_ini[1]: ", 
        ytopo[1] - ytopo_marker_chain_ini[1])

    ytopo_regional = (
        deflections
        .+ ytopo_marker_chain_ini
        .+ (ytopo[1] - ytopo_marker_chain_ini[1])
    )
    return ytopo_regional, xcoor_blocks
end

end # module 