module SimpleSedimentationManager

import Printf: @sprintf
import Statistics: mean
import EarthBox.ModelStructureManager: ModelStructure
import EarthBox.ModelStructureManager: calculate_structure!
import EarthBox.ModelStructureManager: get_basic_structure
import EarthBox.Compaction.CompactionTools: compact_or_decompact
import EarthBox.PlotToolsManager.Charts: plot_ncurves
import EarthBox.PlotToolsManager.Charts: make_plot_name
import EarthBox.EarthBoxDtypes: TestInfoDictType
import EarthBox.EarthBoxDtypes: ModelInputDictType
import EarthBox.EarthBoxDtypes: MaterialsDictType
import ..GetData: get_inputs
import ...TestResults
import ...BenchmarksStruct: Benchmarks

function compare_numerical_to_analytical(
    bench::Benchmarks
)::Tuple{Vector{Union{String, Float64}}, String}
    simple_sedimentation = SimpleSedimentation(bench)
    make_benchmark_plot(simple_sedimentation)
    return get_test_results(simple_sedimentation)
end

mutable struct SimpleSedimentation
    itime_step::Int
    main_paths::Dict{String, String}
    test_info_dict::TestInfoDictType
    model_dict::ModelInputDictType
    materials_dict::Dict{Any, Any}
    materials_library_dict::MaterialsDictType
    tmyr::Float64
    relative_error_limit_percentage::Float64
    gridx_b::Vector{Float64}
    thickness_numerical::Vector{Float64}
    thickness_analytical::Vector{Float64}
end

function SimpleSedimentation(bench::Benchmarks)
    input_dict = get_inputs(bench.itime_step, bench.main_paths)
    
    model_structure = ModelStructure(
        input_dict["marknum"],
        input_dict["marker_matid"],
        input_dict["marker_x_m"],
        input_dict["marker_y_m"],
        input_dict["xnum"],
        input_dict["gridx_b"],
        input_dict["xstp_avg"],
        input_dict["ysize"],
        Int16[1, 2], # matids sticky
        Int16[2], # matids sticky water
        Int16[3, 4], # matids rock
        Int16[-1], # matids crust
        Int16[-1], # matids lithosphere
        Int16[4], # matids asthenosphere
        Int16[3] # matids sediments
    )
    
    calculate_structure!(model_structure)
    
    (
        _ytopo_numerical,
        _ymoho_numerical,
        _ytopo_smooth_numerical,
        _ylith_base_numerical,
        _thick_water_numerical,
        _thick_crust_numerical,
        _thick_mantle_lith_numerical,
        _thick_asthenosphere_numerical,
        thick_sediments_numerical
    ) = get_basic_structure(model_structure)
    
    thickness_analytical = calculate_analytical_thickness()
    thick_sediments_analytical = fill(thickness_analytical, input_dict["xnum"])
    
    return SimpleSedimentation(
        bench.itime_step,
        bench.main_paths,
        bench.test_info_dict,
        bench.model_dict,
        bench.materials_dict,
        bench.materials_library_dict,
        input_dict["tmyr"],
        10.0,
        input_dict["gridx_b"],
        thick_sediments_numerical,
        thick_sediments_analytical
    )
end

function calculate_analytical_thickness()::Float64
    nstep_model = 4
    pelagic_sedimentation_rate_mm_yr = 20.0
    timestep_viscoelastic_yr = 50_000.0
    timestep_transport_yr = 10_000.0
    nstep_transport = floor(Int, timestep_viscoelastic_yr/timestep_transport_yr)
    
    porosity_initial = 0.4
    depth_decay_term = 1.0/2500.0
    
    thickness_sediment_per_timestep_meters = (
        pelagic_sedimentation_rate_mm_yr * timestep_viscoelastic_yr / 1000.0
    )
    
    depositional_thickness = thickness_sediment_per_timestep_meters / nstep_transport
    
    thickness_total = 0.0
    
    for j in 1:nstep_model
        for i in 1:nstep_transport
            if j == 1 && i == 1
                thickness_total = depositional_thickness
            else
                top_initial = 0.0
                bottom_initial = thickness_total
                top_new = depositional_thickness
                thickness_new = compact_or_decompact(
                    porosity_initial,
                    depth_decay_term,
                    top_initial,
                    bottom_initial,
                    top_new
                )
                thickness_total = thickness_new + depositional_thickness
            end
        end
    end
    
    return thickness_total
end

function make_benchmark_plot(data::SimpleSedimentation)::Nothing
    data_xy = Dict{String, Any}(
        "x_arrays" => [
            data.gridx_b ./ 1000.0,
            data.gridx_b ./ 1000.0
        ],
        "y_arrays" => [
            data.thickness_numerical,
            data.thickness_analytical
        ],
        "labels" => [
            "Sediment Thickness (Numerical)",
            "Sediment Thickness (Analytical)"
        ],
        "colors" => [:black, :red],
        "line_styles" => [:solid, :solid],
        "line_widths" => [1, 1],
        "line_colors" => [:black, :red],
        "marker_sizes" => [5, 5],
        "marker_edge_colors" => [:black, :black],
        "marker_edge_widths" => [0.5, 0.5],
        "fill_styles" => [:none, :none]
    )
    
    axis_labels = ["x (km)", "Sediment Thickness (m)"]
    plot_dimensions_xy = [0.0, 140.0, 0.0, 6000.0]
    
    boxtext_x = plot_dimensions_xy[2] * 0.65
    boxtext_y = plot_dimensions_xy[4] * 0.98
    boxtext = ""
    boxtext_info = [boxtext_x, boxtext_y, boxtext]
    
    title_base_text = "Simple Sedimentation with Compaction Benchmark"
    title = "$(title_base_text): Numerical Time $(@sprintf("%.2f", data.tmyr)) Myr"
    
    plot_name = make_plot_name(
        data.main_paths["model_name"] * "_topo",
        data.itime_step
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
        "xtick_size" => 20.0,
        "ytick_size" => 1000.0
    )
    
    plot_ncurves(chart_input)
    return nothing
end

function get_test_results(
    data::SimpleSedimentation
)::Tuple{Vector{Union{String, Float64}}, String}
    numerical_solution = [mean(data.thickness_numerical)]
    analytical_solution = [data.thickness_analytical[1]]
    
    result, result_msg = TestResults.get_test_results_numerical_vs_analytical(
        numerical_solution,
        analytical_solution,
        data.relative_error_limit_percentage / 100.0
    )
    return result, result_msg
end

end # module 