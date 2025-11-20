module PlasticityBenchmarkK10

import Plots
import Plots: gr
import Printf: @sprintf
import EarthBox.InputTools.Reader: make_parameters_dict
import EarthBox.PlotToolsManager.Charts: plot_ncurves
import EarthBox.PlotToolsManager.Charts: make_plot_name
import EarthBox.PlotToolsManager.Charts: check_output_directory
import EarthBox.EarthBoxDtypes: TestInfoDictType
import EarthBox.EarthBoxDtypes: ModelInputDictType
import EarthBox.EarthBoxDtypes: MaterialsDictType
import EarthBox.PlotSettingsManager: PLOT_SETTINGS
import EarthBox.StaggeredGrid.Spacing: initialize_spacing_arrays
import EarthBox.StaggeredGrid.Spacing: update_basic_1Dgrid_spacing!
import EarthBox.StaggeredGrid.VxGridTools: vxgrid_y_or_z_coordinates_loop!
import EarthBox.StaggeredGrid.VxGridTools: vxgrid_y_or_z_spacing_loop!
import EarthBox.StaggeredGrid.VyGridTools: vygrid_x_or_z_coordinates_loop!
import EarthBox.StaggeredGrid.VyGridTools: vygrid_x_or_z_spacing_loop!
import EarthBox.Markers.MarkerMaterials.InitManager.WeakNotch: calculate_seed_limits_for_weak_notch
import EarthBox.Interpolation.GridToMarker: get_marker_value_from_pressure_grid_without_mapping_input
import ...Reader
import ...BenchmarkTools
import ...TestResults
import ...BenchmarksStruct: Benchmarks

function compare_numerical_to_numerical(
    bench::Benchmarks
)::Tuple{Vector{Union{String, Float64}}, String}
    plasticity_benchmark = PlasticityBenchmarkKaus10(bench)
    make_benchmark_plot(plasticity_benchmark)
    return get_test_results(plasticity_benchmark)
end

mutable struct PlasticityBenchmarkKaus10
    itime_step::Int
    main_paths::Dict{String, String}
    test_info_dict::TestInfoDictType
    model_dict::ModelInputDictType
    materials_dict::Dict{Any, Any}
    materials_library_dict::MaterialsDictType
    tmyr::Float64
    relative_error_limit_percentage::Float64
    input_dict::Dict{String, Any}
    shear_zone_dict::Dict{String, Any}
    inclusion_thick::Float64
    inclusion_x::Vector{Float64}
    inclusion_y::Vector{Float64}
    left_fault_stick_x::Vector{Float64}
    left_fault_stick_y::Vector{Float64}
    right_fault_stick_x::Vector{Float64}
    right_fault_stick_y::Vector{Float64}
    fault_dip_array_numerical::Vector{Float64}
    fault_dip_array_target::Vector{Float64}
    fault_dip_array_coulomb::Vector{Float64}
    fault_dip_array_roscoe::Vector{Float64}
    fault_dip_array_arthur::Vector{Float64}
end

function PlasticityBenchmarkKaus10(bench::Benchmarks)
    input_dict = get_inputs(bench)
    shear_zone_dict = calc_x_locations_of_shear_zone_calculation_points(input_dict)
    inclusion_x = zeros(Float64, 5)
    inclusion_y = zeros(Float64, 5)
    inclusion_thick = fill_inclusion_arrays!(shear_zone_dict, inclusion_x, inclusion_y)
    
    # Initialize arrays
    left_fault_stick_x = zeros(Float64, 2)
    left_fault_stick_y = zeros(Float64, 2)
    right_fault_stick_x = zeros(Float64, 2)
    right_fault_stick_y = zeros(Float64, 2)
    fault_dip_array_numerical = zeros(Float64, 1)
    fault_dip_array_target = zeros(Float64, 1)
    fault_dip_array_coulomb = zeros(Float64, 1)
    fault_dip_array_roscoe = zeros(Float64, 1)
    fault_dip_array_arthur = zeros(Float64, 1)
    
    plasticity_benchmark = PlasticityBenchmarkKaus10(
        bench.itime_step, bench.main_paths, bench.test_info_dict, bench.model_dict,
        bench.materials_dict, bench.materials_library_dict, input_dict["tmyr"], 4.0,
        input_dict, shear_zone_dict, inclusion_thick, inclusion_x, inclusion_y,
        left_fault_stick_x, left_fault_stick_y, right_fault_stick_x, right_fault_stick_y,
        fault_dip_array_numerical, fault_dip_array_target, fault_dip_array_coulomb,
        fault_dip_array_roscoe, fault_dip_array_arthur
    )
    
    calculate_fault_geometry(plasticity_benchmark)
    return plasticity_benchmark
end

function get_inputs(bench::Benchmarks)::Dict{String, Any}
    file_base_name = "Sxy_MPa"
    gridx_b_km, gridy_b_km, sxy_mpa, tmyr = Reader.read_scalar_grid(
        bench.main_paths["post_proc_input_path"],
        file_base_name, bench.itime_step
        )
    
    sxy = sxy_mpa * 1e6
    gridx_b = gridx_b_km * 1000.0
    gridy_b = gridy_b_km * 1000.0
    xnum = length(gridx_b_km)
    ynum = length(gridy_b_km)
    xsize = gridx_b[xnum]
    ysize = gridy_b[ynum]
    
    xstp_b, ystp_b = initialize_spacing_arrays(gridx_b, gridy_b)
    update_basic_1Dgrid_spacing!(xnum, gridx_b, xstp_b)
    update_basic_1Dgrid_spacing!(ynum, gridy_b, ystp_b)
    
    gridy_vx = zeros(Float64, ynum + 1)
    ystp_vx = zeros(Float64, ynum)
    gridx_vy = zeros(Float64, xnum + 1)
    xstp_vy = zeros(Float64, xnum)
    
    vxgrid_y_or_z_coordinates_loop!(ynum, gridy_b, ystp_b, gridy_vx)
    vxgrid_y_or_z_spacing_loop!(ynum, gridy_b, ystp_b, ystp_vx)
    vygrid_x_or_z_coordinates_loop!(xnum, gridx_b, xstp_b, gridx_vy)
    vygrid_x_or_z_spacing_loop!(xnum, gridx_b, xstp_b, xstp_vy)
    
    file_base_name = "Eii_log"
    pr_gridx_km, pr_gridy_km, eii_log, tmyr = Reader.read_scalar_grid(
        bench.main_paths["post_proc_input_path"],
        file_base_name, bench.itime_step
        )
    
    eii = 10.0 .^ eii_log
    
    parameters_dict = make_parameters_dict(bench.model_dict)
    thick_air = parameters_dict["thick_air"][1]
    w_seed = parameters_dict["w_seed"][1]
    x_seed = parameters_dict["x_seed"][1]
    
    matid = 3
    material_name = bench.materials_dict[matid]["mat_name"][1]
    material_dict_layer1 = bench.materials_library_dict[material_name]
    friction_angle = material_dict_layer1["friction_angle_initial"]
    
    dilation_angle = 0.0
    
    coulomb_angle = 45.0 + friction_angle / 2.0
    roscoe_angle = 45.0 + dilation_angle / 2.0
    arthur_angle = 45.0 + (friction_angle + dilation_angle) / 4.0
    
    return Dict{String, Any}(
        "tmyr" => tmyr,
        "gridx_b" => gridx_b,
        "gridy_b" => gridy_b,
        "xstp_b" => xstp_b,
        "ystp_b" => ystp_b,
        "pr_gridx" => pr_gridx_km * 1000.0,
        "pr_gridy" => pr_gridy_km * 1000.0,
        "gridy_vx" => gridy_vx,
        "gridx_vy" => gridx_vy,
        "ystp_vx" => ystp_vx,
        "xstp_vy" => xstp_vy,
        "xnum" => xnum,
        "ynum" => ynum,
        "xsize" => xsize,
        "ysize" => ysize,
        "thick_air" => thick_air,
        "sxy" => sxy,
        "eii" => eii,
        "eii_log" => eii_log,
        "friction_angle" => friction_angle,
        "w_seed" => w_seed,
        "x_seed" => x_seed,
        "coulomb_angle" => coulomb_angle,
        "roscoe_angle" => roscoe_angle,
        "arthur_angle" => arthur_angle
    )
end

function fill_inclusion_arrays!(
    shear_zone_dict::Dict{String, Any},
    inclusion_x::Vector{Float64},
    inclusion_y::Vector{Float64}
)::Float64
    xmin_inclusion = shear_zone_dict["xmin_inclusion"]
    xmax_inclusion = shear_zone_dict["xmax_inclusion"]
    ymin_inclusion = shear_zone_dict["ymin_inclusion"]
    ymax_inclusion = shear_zone_dict["ymax_inclusion"]
    
    println("xmin_inclusion : ", xmin_inclusion)
    println("xmax_inclusion : ", xmax_inclusion)
    println("ymin_inclusion : ", xmax_inclusion)
    println("ymax_inclusion : ", ymax_inclusion)
    
    inclusion_x[1] = xmin_inclusion
    inclusion_y[1] = ymin_inclusion
    inclusion_x[2] = xmax_inclusion
    inclusion_y[2] = ymin_inclusion
    inclusion_x[3] = xmax_inclusion
    inclusion_y[3] = ymax_inclusion
    inclusion_x[4] = xmin_inclusion
    inclusion_y[4] = ymax_inclusion
    inclusion_x[5] = xmin_inclusion
    inclusion_y[5] = ymin_inclusion
    
    return ymax_inclusion - ymin_inclusion
end

function calculate_fault_geometry(data::PlasticityBenchmarkKaus10)
    x1_left = data.shear_zone_dict["x1_left"]
    x2_left = data.shear_zone_dict["x2_left"]
    x1_right = data.shear_zone_dict["x1_right"]
    x2_right = data.shear_zone_dict["x2_right"]
    
    y1_left, y2_left, y1_right, y2_right = calculate_y_locations_of_strain_rate_maxima(
        data.input_dict, x1_left, x2_left, x1_right, x2_right, data.inclusion_thick
    )
    
    data.left_fault_stick_x[1] = x1_left
    data.left_fault_stick_x[2] = x2_left
    data.left_fault_stick_y[1] = y1_left
    data.left_fault_stick_y[2] = y2_left
    data.right_fault_stick_x[1] = x1_right
    data.right_fault_stick_x[2] = x2_right
    data.right_fault_stick_y[1] = y1_right
    data.right_fault_stick_y[2] = y2_right
    
    dy = abs(y2_left - y1_left)
    dx = abs(x2_left - x1_left)
    theta_left = atan(dy/dx) * 180.0 / π
    dy = abs(y2_right - y1_right)
    dx = abs(x2_right - x1_right)
    theta_right = atan(dy/dx) * 180.0 / π
    theta = (theta_left + theta_right) / 2.0
    data.fault_dip_array_numerical[1] = theta
    
    data.fault_dip_array_coulomb[1] = data.input_dict["coulomb_angle"]
    data.fault_dip_array_roscoe[1] = data.input_dict["roscoe_angle"]
    data.fault_dip_array_arthur[1] = data.input_dict["arthur_angle"]
    
    data.fault_dip_array_target[1] = (
        (data.input_dict["arthur_angle"] + data.input_dict["coulomb_angle"]) / 2.0
    )
    
    println("")
    println("theta_left: ", theta_left)
    println("theta_right: ", theta_right)
    println("theta: ", theta)
    println("coulomb_angle: ", data.input_dict["coulomb_angle"])
    println("roscoe_angle: ", data.input_dict["roscoe_angle"])
    println("arthur_angle: ", data.input_dict["arthur_angle"])
    println("")
end

function make_benchmark_plot(data::PlasticityBenchmarkKaus10)::Nothing
    # Turn off interactive plotting
    gr()
    # Create figure and axes
    p = Plots.plot(
        xlims=(0, data.input_dict["xsize"]),
        ylims=(0.0, data.input_dict["ysize"]),
        xticks=[0, 10000, 20000, 30000, 40000],
        yticks=[0, 2000, 4000, 6000, 8000, 10000, 12000],
        yflip=true,
        aspect_ratio=:equal,
        legendfontsize=8,
        titlefontsize=8,
        guidefontsize=8,
        tickfontsize=6,
        size=(1200, 360) # Equivalent to figsize=(40,12)
    )
    # Create heatmap
    Plots.heatmap!(
        p,
        data.input_dict["pr_gridx"],
        data.input_dict["pr_gridy"], 
        data.input_dict["eii_log"],
        c=:bwr,
        clims=(-18, -13),
        colorbar=:horizontal,
        colorbar_ticks=[-18, -17, -16, -15, -14, -13],
        colorbar_title="",
        colorbar_titlefontsize=8,
        colorbar_tickfontsize=6
    )
    # Plot fault lines and inclusion
    Plots.plot!(
        p,
        data.left_fault_stick_x,
        data.left_fault_stick_y,
        label="Left Shear Zone",
        color=:black,
        linestyle=:solid,
        marker=:circle,
        markersize=4,
        markerstrokecolor=:black,
        markerstrokewidth=0.5,
        markerfillcolor=nothing
    )
    
    Plots.plot!(
        p,
        data.right_fault_stick_x,
        data.right_fault_stick_y,
        label="Right Shear Zone", 
        color=:black,
        linestyle=:solid,
        marker=:circle,
        markersize=4,
        markerstrokecolor=:black,
        markerstrokewidth=0.5,
        markerfillcolor=nothing
    )
    Plots.plot!(
        p,
        data.inclusion_x,
        data.inclusion_y,
        label="Inclusion",
        color=:black,
        linestyle=:solid,
        marker=:circle,
        markersize=1,
        markerstrokecolor=:black,
        markerstrokewidth=0.5,
        markerfillcolor=nothing
    )
    
    # Set title and labels
    Plots.title!("Plasticity Benchmark Kaus (2010)")
    Plots.xlabel!("X (m)")
    Plots.ylabel!("Y (m)")
    
    plot_name = "plasticity_benchmark_kaus10.$(PLOT_SETTINGS.plot_extension)"
    plot_file_path = joinpath(data.main_paths["post_proc_output_path"], plot_name)
    println("plot_file_path : ", plot_file_path)
    check_output_directory(dirname(plot_file_path))
    Plots.savefig(p, plot_file_path)
    return nothing
end

function get_test_results(
    data::PlasticityBenchmarkKaus10
)::Tuple{Vector{Union{String, Float64}}, String}
    result, result_msg = TestResults.get_test_results_numerical_vs_analytical(
        data.fault_dip_array_numerical,
        data.fault_dip_array_target,
        data.relative_error_limit_percentage / 100.0
    )
    return result, result_msg
end

function calc_x_locations_of_shear_zone_calculation_points(
    input_dict::Dict{String, Any}
)::Dict{String, Any}
    xmin_inclusion, xmax_inclusion, ymin_inclusion, ymax_inclusion = 
        calculate_seed_limits_for_weak_notch(
            input_dict["x_seed"],
            input_dict["ysize"],
            input_dict["w_seed"]
        )
    
    ysize = input_dict["ysize"]
    println("ysize : ", ysize)
    println("xmin_inclusion : ", xmin_inclusion)
    println("xmax_inclusion : ", xmax_inclusion)
    println("ymin_inclusion : ", xmax_inclusion)
    println("ymax_inclusion : ", ymax_inclusion)
    
    x1_left = xmin_inclusion - input_dict["w_seed"] / 4.0
    x2_left = xmin_inclusion - (input_dict["w_seed"] / 4.0 + 2000.0)
    x1_right = xmax_inclusion + input_dict["w_seed"] / 4.0
    x2_right = xmax_inclusion + (input_dict["w_seed"] / 4.0 + 2000.0)
    
    return Dict{String, Any}(
        "x1_left" => x1_left,
        "x2_left" => x2_left,
        "x1_right" => x1_right,
        "x2_right" => x2_right,
        "xmin_inclusion" => xmin_inclusion,
        "xmax_inclusion" => xmax_inclusion,
        "ymin_inclusion" => ymin_inclusion,
        "ymax_inclusion" => ymax_inclusion
    )
end

function calculate_y_locations_of_strain_rate_maxima(
    input_dict::Dict{String, Any},
    x1_left::Float64,
    x2_left::Float64,
    x1_right::Float64,
    x2_right::Float64,
    inclusion_thick::Float64
)::Tuple{Float64, Float64, Float64, Float64}
    ysize = input_dict["ysize"]
    thick_air = input_dict["thick_air"]
    y_initial = thick_air
    y_thick = ysize - y_initial - inclusion_thick
    
    ynum = input_dict["ynum"]
    npoints = ynum * 10
    delta_y = y_thick / (npoints - 1)
    
    y1_left, _ = find_local_maximum_in_y_direction(
        npoints, x1_left, delta_y, y_initial,
        input_dict["gridx_b"], input_dict["gridy_b"],
        input_dict["xstp_b"], input_dict["ystp_b"],
        input_dict["eii"], input_dict["gridy_vx"],
        input_dict["ystp_vx"], input_dict["gridx_vy"],
        input_dict["xstp_vy"]
    )
    
    y2_left, _ = find_local_maximum_in_y_direction(
        npoints, x2_left, delta_y, y_initial,
        input_dict["gridx_b"], input_dict["gridy_b"],
        input_dict["xstp_b"], input_dict["ystp_b"],
        input_dict["eii"], input_dict["gridy_vx"],
        input_dict["ystp_vx"], input_dict["gridx_vy"],
        input_dict["xstp_vy"]
    )
    
    y1_right, _ = find_local_maximum_in_y_direction(
        npoints, x1_right, delta_y, y_initial,
        input_dict["gridx_b"], input_dict["gridy_b"],
        input_dict["xstp_b"], input_dict["ystp_b"],
        input_dict["eii"], input_dict["gridy_vx"],
        input_dict["ystp_vx"], input_dict["gridx_vy"],
        input_dict["xstp_vy"]
    )
    
    y2_right, _ = find_local_maximum_in_y_direction(
        npoints, x2_right, delta_y, y_initial,
        input_dict["gridx_b"], input_dict["gridy_b"],
        input_dict["xstp_b"], input_dict["ystp_b"],
        input_dict["eii"], input_dict["gridy_vx"],
        input_dict["ystp_vx"], input_dict["gridx_vy"],
        input_dict["xstp_vy"]
    )
    
    return y1_left, y2_left, y1_right, y2_right
end

function find_local_maximum_in_y_direction(
    npoints::Int,
    x_point::Float64,
    delta_y::Float64,
    y_initial::Float64,
    gridx_b::Vector{Float64},
    gridy_b::Vector{Float64},
    xstp_b::Vector{Float64},
    ystp_b::Vector{Float64},
    scalar_grid_array::Matrix{Float64},
    gridy_vx::Vector{Float64},
    ystp_vx::Vector{Float64},
    gridx_vy::Vector{Float64},
    xstp_vy::Vector{Float64}
)::Tuple{Float64, Float64}
    y_max = -1e32
    scalar_max = -1e32
    
    for i in 1:npoints
        y_point = y_initial + delta_y * (i - 1)
        scalar_value = get_marker_value_from_pressure_grid_without_mapping_input(
            gridx_b, gridy_b, xstp_b, ystp_b,
            gridy_vx, ystp_vx, gridx_vy, xstp_vy,
            scalar_grid_array, x_point, y_point
        )
        if scalar_value > scalar_max
            scalar_max = scalar_value
            y_max = y_point
        end
    end
    
    return y_max, scalar_max
end

end # module 