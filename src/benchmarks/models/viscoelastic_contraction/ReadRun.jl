""" Run model by reading input yaml files.
"""
module ReadRun

import EarthBox
import EarthBox.PlotToolsManager: ModelPlots2DManager
import EarthBox.PlotToolsManager.ModelPlots2DManager: GridPlotsManager
import EarthBox.PlotToolsManager.ModelPlots2DManager: MarkerPlotsManager
import EarthBox.PlotToolsManager.ModelPlots2DManager: RheologyPlotsManager
import EarthBox.MaterialLibraryCollection: MaterialLibrary

function get_model_path()::String
    return joinpath(@__DIR__, "model.yml")
end

function get_materials_path()::String
    return joinpath(@__DIR__, "materials.yml")
end

function get_material_lib_path()::String
    return MaterialLibrary().benchmarks.path
end

function get_output_path()::String
    return "/mnt/extradrive1/earthbox_benchmark_results_2025-10-07_10-08-47/viscoelastic_contraction_output"
end

function run_model(;
    output_dir::String = get_output_path(),
    make_backup::Bool = false,
    restart_from_backup::Bool = false,
    use_mumps::Bool = false,
    nproc::Int = 1
)::Nothing
    eb = EarthBox.EarthBoxState(
        restart_from_backup = restart_from_backup,
        paths = Dict(
            "model_input_file" => get_model_path(),
            "materials_input_file" => get_materials_path(),
            "materials_library_file" => get_material_lib_path(),
            "output_dir" => output_dir
        )
    )
    # Activate marker age output
    eb.model_manager.config.output.marker_output["marker_age"] = true

    eb.model_manager.config.solver.use_mumps = use_mumps
    eb.model_manager.config.solver.nproc = nproc
    EarthBox.ModelManager.initialize_model!(eb.model_manager)
    EarthBox.run_time_steps(eb, make_backup=make_backup)
    return nothing
end

function get_model_plots_2d(;output_dir::String = get_output_path())::ModelPlots2DManager.ModelPlots2D
    return ModelPlots2DManager.ModelPlots2D(
        plot_output_path=joinpath(output_dir, "plots"),
        material_library_file_path=get_material_lib_path(),
        material_model_file_path=get_materials_path(),
        model_output_path=output_dir,
        nsteps=50,
        nskip=0,
        istart=0,
        iplot_contour_labels=1,
        contour_label_fontsize=5,
        color_map="bwr",
        figure_dpi=150.0,
        figsize=(10.0, 5.0),
        length_units="km",
        time_units="Myr",
        velocity_units="cm/yr",
        dimensions=(0.0, 1000.0, 0.0, 160.0),
        use_close_up=false,
        dim_close_up=(50.0, 100.0, 50.0, 100.0),
        xyspacing=(1.0, 20.0),
        linewidth=0.01,
        title_fontsize=12,
        legend_fontsize=8,
        axis_title_fontsize=8,
        axis_labels_fontsize=8,
        axis_ticks_fontsize=8
    )
end

function plot_yield_stress(;output_dir::String = get_output_path())::Nothing
    mp2d = get_model_plots_2d(output_dir=output_dir)

    RheologyPlotsManager.plot_yield_strength(
        mp2d.rheology_plots,
        depth_plot_limit=150.0,
        depth_axis_spacing=10.0,
        temperature_plot_limit=1400.0,
        temperature_axis_spacing=200.0,
        yield_stress_plot_limit=600.0,
        yield_stress_axis_spacing=50.0,
        strain_rate=1e-15,
        thickness_upr_cont_crust_meters=22_000.0,
        thickness_lwr_cont_crust_meters=13_000.0,
        thickness_lithosphere_meters=125_000.0,
        thickness_asthenosphere_meters=25_000.0,
        dy_meters=100.0,
        expansivity=3e-5,
        compressibility=1e-11,
        density_upper_continental_crust=2822.72,
        density_lower_continental_crust=2900.0,
        density_mantle_lithosphere=3250.0,
        density_asthenosphere=3300.0,
        iuse_linear_segments=false,
        temperature_top_celsius=0.0,
        temperature_moho_celsius=600.0,
        temperature_base_lith_celsius=1330.0,
        adiabatic_gradient_kelvin_km=0.4,
        conductivity_upper_crust=2.25,
        conductivity_lower_crust=2.2,
        conductivity_mantle=2.0,
        heat_production_upper_crust=1.8e-6,
        heat_production_lower_crust=0.5e-6,
        heat_production_mantle=0.0,
        thickness_thermal_lithosphere=125_000.0,
        iuse_fluid_pressure_for_yield=1,
        log10_viscosity_plot_min=20.0,
        log10_viscosity_plot_max=26.0,
        extension=".pdf",
        figsize=(3.0, 6.0)
        )
    return nothing
end

function plot_markers(;output_dir::String = get_output_path())::Nothing
    println("Plotting markers")
    mp2d = get_model_plots_2d(output_dir=output_dir)
    
    MarkerPlotsManager.plot_markers(
        mp2d.marker_plots,
        plot_type="Composition",#"CompositionGravity", #"CompositionHeatFlow",
        marker_size=2.0,
        decimation_factor=4,
        plot_mesh=1,
        mesh_line_width=0.1,
        plot_contour_labels=0,
        contour_line_width=0.2,
        contour_line_color="black",
        plot_temperature_contours=0,
        temperature_min=0.0001,
        temperature_max=1380.0,
        temperature_contour_interval=100.0,
        temperature_contour_color="black",
        plot_topography=1,
        topo_line_width=1.0,
        topo_line_color="red",
        plot_base_level=1,
        base_level_line_width=1.0,
        base_level_line_color="blue",
        plot_plastic_strain=false,
        strain_min=1.0,
        strain_max=6.0,
        strain_contour_interval=0.25,
        strain_cmap="inferno",
        heatflow_min=0.0,
        heatflow_max=200.0,
        heatflow_spacing=50.0,
        plot_intrusive_age=0,
        age_min_intrusive=0.0,
        age_max_intrusive=10.0,
        age_contour_interval_intrusive=0.2,
        age_cmap_intrusive="rainbow"
    )

end

function plot_scalars()
    mp2d = get_model_plots_2d()
    scalar_names = GridPlotsManager.ScalarNames()

    GridPlotsManager.plot_scalars(
        mp2d.grid_plots,
        scalar_names.viscosity,
        contour_interval=0.5,
        minimum_value=18.0,
        maximum_value=26.0,
        plot_contours=true,
        grid_plot_type="nomesh"
    )

    # GridPlotsManager.plot_scalars(
    #     mp2d.grid_plots,
    #     scalar_names.temperature,
    #     contour_interval=100.0,
    #     minimum_value=0.0,
    #     maximum_value=1300.0,
    #     plot_contours=true,
    #     grid_plot_type="nomesh"
    # )
    
    # GridPlotsManager.plot_scalars(
    #     mp2d.grid_plots,
    #     scalar_names.velocity_x,
    #     contour_interval=0.25,
    #     minimum_value=-2.5,
    #     maximum_value=2.5,
    #     plot_contours=false,
    #     grid_plot_type="nomesh"
    # )
    
    # GridPlotsManager.plot_scalars(
    #     mp2d.grid_plots,
    #     scalar_names.velocity_y,
    #     contour_interval=0.25,
    #     minimum_value=-2.5,
    #     maximum_value=2.5,
    #     plot_contours=false,
    #     grid_plot_type="nomesh"
    # )
    
    # GridPlotsManager.plot_scalars(
    #     mp2d.grid_plots,
    #     scalar_names.velocity_mag,
    #     contour_interval=0.25,
    #     minimum_value=-2.5,
    #     maximum_value=2.5,
    #     plot_contours=false,
    #     grid_plot_type="nomesh"
    # )
    
end

# function plot_velocity()
#     mp2d = get_model_plots_2d()
#     mp2d.plot_velocity()
# end
# 
# function plot_markers()
#     mp2d = get_model_plots_2d()
#     plot_types = eb.model_manager.model_plots_2d.marker_plot_names
#     mp2d.plot_markers(
#         plot_type=plot_types.composition,
#         plot_topography=true,
#         plot_mesh=false,
#         marker_size=0.25
#     )
# end
# 
# """Get function based on option name."""
# function get_func(option_name::String)::Function
#     func_dict = Dict(
#         "run_model" => run_model,
#         "plot_markers" => plot_markers,
#         "plot_scalars" => plot_scalars,
#         "plot_velocity" => plot_velocity
#     )
#     
#     if !haskey(func_dict, option_name)
#         throw(ArgumentError(
#             "$(option_name) is not a valid option. " *
#             "Valid options are $(keys(func_dict))"
#         ))
#     end
#     return func_dict[option_name]
# end

end # module