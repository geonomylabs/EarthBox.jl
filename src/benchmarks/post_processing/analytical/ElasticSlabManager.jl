module ElasticSlabManager

import Printf
import Plots
import EarthBox.InputTools.Reader: make_parameters_dict
import EarthBox.PlotToolsManager.Charts: plot_ncurves
import EarthBox.PlotToolsManager.Charts: make_plot_name
import EarthBox.PlotToolsManager.Charts: check_output_directory
import EarthBox.EarthBoxDtypes: TestInfoDictType
import EarthBox.EarthBoxDtypes: ModelInputDictType
import EarthBox.EarthBoxDtypes: MaterialsDictType
import EarthBox.Markers.MarkerMaterials.InitManager.ElasticSlab: in_elastic_slab
import EarthBox.JLDTools: get_marker_coordinate_arrays
import ...BenchmarkTools
import ...TestResults
import ...BenchmarksStruct: Benchmarks

function compare_elastic_slab_marker_position_to_initial(
    bench::Benchmarks
)::Tuple{Vector{Union{String, Float64}}, String}
    elastic_slab = ElasticSlab(bench)
    make_benchmark_plot(elastic_slab)
    return get_test_results(elastic_slab)
end

mutable struct ElasticSlab
    itime_step::Int
    main_paths::Dict{String, String}
    test_info_dict::TestInfoDictType
    model_dict::ModelInputDictType
    materials_dict::Dict{Any, Any}
    materials_library_dict::MaterialsDictType
    tmyr::Float64
    relative_error_limit_percentage::Float64
    markers_x_initial_m::Vector{Float64}
    markers_y_initial_m::Vector{Float64}
    markers_x_final_m::Vector{Float64}
    markers_y_final_m::Vector{Float64}
    markers_x_slab_m::Vector{Float64}
    markers_y_slab_m::Vector{Float64}
    markers_ids_slab::Vector{Int64}
    marker_dist_from_origin_initial::Vector{Float64}
    marker_dist_from_origin_final::Vector{Float64}
    difference::Vector{Float64}
end

function ElasticSlab(bench::Benchmarks)
    _, markers_x_initial_m, markers_y_initial_m = get_marker_coordinate_arrays(
        1, bench.main_paths["post_proc_input_path"]
    )

    tmyr, markers_x_final_m, markers_y_final_m = get_marker_coordinate_arrays(
        bench.itime_step, bench.main_paths["post_proc_input_path"]
    )

    parameters_dict = make_parameters_dict(bench.model_dict)
    xsize = parameters_dict["xsize"][1]
    ysize = parameters_dict["ysize"][1]

    markers_ids_slab, markers_x_slab_m, markers_y_slab_m = 
        get_marker_ids_of_elastic_slab(
            markers_x_initial_m, markers_y_initial_m, xsize, ysize
        )

    marker_dist_from_origin_initial = calc_marker_dist_from_origin(
        markers_x_initial_m, markers_y_initial_m, markers_ids_slab
    )
    marker_dist_from_origin_final = calc_marker_dist_from_origin(
        markers_x_final_m, markers_y_final_m, markers_ids_slab
    )
    difference = abs.(
        marker_dist_from_origin_final .- marker_dist_from_origin_initial
    )

    return ElasticSlab(
        bench.itime_step, bench.main_paths, bench.test_info_dict, bench.model_dict,
        bench.materials_dict, bench.materials_library_dict, tmyr, 0.3,
        markers_x_initial_m, markers_y_initial_m, markers_x_final_m, markers_y_final_m,
        markers_x_slab_m, markers_y_slab_m, markers_ids_slab,
        marker_dist_from_origin_initial, marker_dist_from_origin_final, difference
    )
end

function make_benchmark_plot(data::ElasticSlab)::Nothing
    figsize = (10, 10)
    figsize_pixels = (figsize[1] * 100, figsize[2] * 100)
   
    axis_labels = ["X Axis (km)", "Y Axis (km)"]

    scatter_plot_obj = Plots.scatter(
        data.markers_x_slab_m ./ 1000.0,
        data.markers_y_slab_m ./ 1000.0,
        marker_z=data.difference,
        clims=(0.0, 2000.0),
        color=:viridis,
        xlabel=axis_labels[1],
        ylabel=axis_labels[2],
        size=figsize_pixels,
        dpi=150,
        markershape=:rect,
        markersize=4.0,
        aspect_ratio=1.0,
        markerstrokewidth=0.0,
        markerstrokecolor=:transparent,
        legend=false
    )
    Plots.plot!(
        scatter_plot_obj, colorbar=true, 
        colorbar_title="Difference between initial and final location (meters)"
        )
    Plots.yflip!(scatter_plot_obj)

    plot_name = "elastic_slab_distance_from_initial_$(data.itime_step).png"
    plot_file_path = joinpath(
        data.main_paths["post_proc_output_path"], plot_name)
    title = string(
        "Marker Distance From Initial Position: ",
        "Time $(Printf.@sprintf("%.2f", data.tmyr)) Myr"
        )
    Plots.title!(scatter_plot_obj, title)
    check_output_directory(dirname(plot_file_path))
    Plots.savefig(scatter_plot_obj, plot_file_path)
    return nothing
end

function get_test_results(
    data::ElasticSlab
)::Tuple{Vector{Union{String, Float64}}, String}
    result, result_msg = TestResults.get_test_results_numerical_vs_analytical(
        data.marker_dist_from_origin_final,
        data.marker_dist_from_origin_initial,
        data.relative_error_limit_percentage / 100.0
    )
    return result, result_msg
end

""" Calculate initial coordinate of reference marker.
"""
function get_marker_ids_of_elastic_slab(
    marker_x_initial::Vector{Float64},
    marker_y_initial::Vector{Float64},
    xsize::Float64,
    ysize::Float64
)::Tuple{Vector{Int64}, Vector{Float64}, Vector{Float64}}
    nmarkers = length(marker_x_initial)
    marker_ids_tmp = fill(-1, nmarkers)
    icount = 0
    
    for imarker in 1:nmarkers
        x_marker = marker_x_initial[imarker]
        y_marker = marker_y_initial[imarker]
        if in_elastic_slab(y_marker, x_marker, ysize, xsize)
            icount += 1
            marker_ids_tmp[icount] = imarker
        end
    end

    marker_ids_slab = zeros(Int64, icount)
    icount = 0
    for i in 1:nmarkers
        test_id = marker_ids_tmp[i]
        if test_id != -1
            icount += 1
            marker_ids_slab[icount] = test_id
        end
    end

    marker_x_slab = zeros(Float64, icount)
    marker_y_slab = zeros(Float64, icount)
    for imarker in 1:icount
        marker_id = marker_ids_slab[imarker]
        marker_x_slab[imarker] = marker_x_initial[marker_id]
        marker_y_slab[imarker] = marker_y_initial[marker_id]
    end

    return marker_ids_slab, marker_x_slab, marker_y_slab
end

""" Calculate the distance of target markers from origin.
"""
function calc_marker_dist_from_origin(
    markers_x::Vector{Float64},
    markers_y::Vector{Float64},
    markers_ids_target::Vector{Int64}
)::Vector{Float64}
    ntarget = length(markers_ids_target)
    marker_dist_from_origin = zeros(Float64, ntarget)
    nmarkers = length(markers_x)
    icount = 0
    
    for imarker in 1:nmarkers
        if imarker in markers_ids_target
            x_marker = markers_x[imarker]
            y_marker = markers_y[imarker]
            dist = sqrt(x_marker^2.0 + y_marker^2.0)
            icount += 1
            marker_dist_from_origin[icount] = dist
        end
    end
    return marker_dist_from_origin
end

end # module 