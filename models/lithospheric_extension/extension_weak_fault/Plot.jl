"""
    Plot.jl

Customized plotting module with user defined plotting functions and command-line
plotter that plot output files stored in the model output directory. Note that 
model output directories are named `<case_name>_output` where `<case_name>` is 
the name of the model case like `case0`, `case1`, `case2`, etc (default=`case0`).

For details on command-line arguments and examples, see the section 
called "Command-line Plotter" or search for `run_cl_plotter` in the EarthBox 
documentation.

# Quick Start

## Plotting via command-line:

```bash
julia Plot.jl <plot_option_name> case_name=<case_name> istart=1 iend=100
```

where `istart` and `iend` are the starting and ending time steps, `<case_name>` 
is the name of the model case and `<plot_option_name>` is the name of the plot 
option name that can be one of the following:

- `marker_plots`: Plot markers.
- `scalar_plots`: Plot scalars.
- `velocity_plots`: Plot velocity.
- `stokes_convergence_plots`: Plot Stokes convergence.
- `yield_strength_plots`: Plot yield strength.

## Plotting via REPL:

```julia
include("Plot.jl")
Plot.marker_plots(
    model_output_path = "/path/to/model_output",
    istart = 1,
    iend = 100
);

Note that `istart` and `iend` are not used by the stokes convergence and yield 
strength plotting functions.

"""
module Plot

using EarthBox

include("Model.jl")
include("Materials.jl")
import .Model: ROOT_PATH_OUTPUT
import .Materials: get_materials_input_dict, MATERIAL_COLLECTION

const dimensions = (0.0, 500.0, 0.0, 160.0)
const xyspacing = (50.0, 10.0)
const model_figsize = (10.0, 5.0)

function scalar_plots(;
    model_output_path::String,
    istart::Int64 = 1,
    iend::Union{Int64, Nothing} = nothing
)::Nothing
    plot_scalars(
        scalar_name = :temperature,
        model_output_path = model_output_path,
        dimensions = dimensions,
        xyspacing = xyspacing,
        istart = istart,
        iend = iend,
        figsize = model_figsize
    )
    plot_scalars(
        scalar_name = :viscosity,
        model_output_path = model_output_path,
        dimensions = dimensions,
        xyspacing = xyspacing,
        istart = istart,
        iend = iend,
        figsize = model_figsize
    )
    return nothing
end

function velocity_plots(;
    model_output_path::String,
    istart::Int64 = 1,
    iend::Union{Int64, Nothing} = nothing
)::Nothing
    plot_velocity(
        model_output_path = model_output_path,
        dimensions        = dimensions,
        xyspacing         = xyspacing,
        istart            = istart,
        iend              = iend,
        decimation_factor = 8,
        scale_factor      = 10.0,
        figsize = model_figsize
    )
    return nothing
end

function marker_plots(;
    model_output_path::String,
    istart::Int64 = 1,
    iend::Union{Int64, Nothing} = nothing
)::Nothing
    plot_markers(
        plot_type=:CompositionHeatFlowGravity,
        model_output_path=model_output_path,
        material_library_file_path=MATERIAL_COLLECTION.path,
        materials_input_dict=get_materials_input_dict(),
        # General plotting parameters
        dimensions=dimensions,
        xyspacing=xyspacing,
        istart=istart, iend=iend,
        figsize = model_figsize,
        xy_location_contour_legend = (1.0, 5.0),
        text_box_font_size = 5,
        # General marker parameters
        # (reducing the decimation factor reduces pixelation of the plot)
        marker_size=2.0, decimation_factor=5,
        plot_mesh=0, mesh_line_width=0.1,
        plot_contour_labels=1, contour_line_width=1.0, contour_line_color=:black,
        # Topography parameters
        plot_topography=1, topo_line_width=1.0, 
        topo_line_color=:red,
        # Base level parameters
        plot_base_level=1, base_level_line_width=1.0,
        base_level_line_color=:blue,
        # Plastic strain parameters
        plot_plastic_strain = true,
        strain_min = 1.0, strain_max = 6.0,
        strain_contour_interval = 0.25, strain_cmap=:inferno,
        # Plastic strain rate parameters (log10(1/s))
        plot_plastic_strain_rate = 1,
        strain_rate_min = -18, strain_rate_max = -12,
        strain_rate_contour_interval = 0.5, strain_rate_cmap = :Reds,
        # Temperature contours parameters (Celsius)
        plot_temperature_contours=1,
        temperature_min=100.0, temperature_max=1400.0,
        temperature_contour_interval=100.0,
        temperature_number_format="%6.0f",
        temperature_contour_color=:black,
        # Heat flow parameters (mW/m^2)
        heatflow_min=0.0, heatflow_max=200.0,
        heatflow_spacing=50.0,
        # Gravity Parameters (mgal)
        gravity_min=-1000.0, gravity_max=1000.0,
        gravity_spacing=200.0
    )
    return nothing
end

function yield_strength_plots(
    ;model_output_path::String,
)::Nothing
    plot_yield_strength(
        model_output_path=model_output_path,
        material_library_file_path=MATERIAL_COLLECTION.path,
        materials_input_dict=get_materials_input_dict(),
        depth_plot_limit=150.0,
        depth_axis_spacing=10.0,
        temperature_plot_limit=1400.0,
        temperature_axis_spacing=200.0,
        yield_stress_plot_limit=600.0,
        yield_stress_axis_spacing=50.0,
        strain_rate=1e-15,
        thickness_upr_cont_crust_meters=22_000.0,
        thickness_lwr_cont_crust_meters=10_000.0,
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

function stokes_convergence_plots(
    ;model_output_path::String,
)::Nothing
    plot_stokes_convergence(
        model_output_path=model_output_path,
        figsize=(16, 4),
        xmin=0,
        xmax=200,
        xspacing=50,
        log_l2_ymin=-7,
        log_l2_ymax=1,
        log_l2_yspacing=1,
        plot_title="TestCase",
        axis_label_size=8,
        axis_title_size=8,
        legend_font_size=8,
        xtick_label_size=6,
        ytick_label_size=6,
        annotation_font_size=4
        )
    return nothing
end

function main()::Nothing
    run_cl_plotter(
        root_path_output = ROOT_PATH_OUTPUT,
        marker_plots_func = marker_plots,
        scalar_plots_func = scalar_plots,
        velocity_plots_func = velocity_plots,
        stokes_convergence_plots_func = stokes_convergence_plots,
        yield_strength_plots_func = yield_strength_plots,
    )
    return nothing
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    Plot.main()
end