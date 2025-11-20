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

const dimensions = (150.0, 350.0, 0.0, 50.0)
const xyspacing = (50.0, 10.0)
const model_figsize = (10.0, 3.0)

function scalar_plots(;
    model_output_path::String,
    istart::Int64 = 1,
    iend::Union{Int64, Nothing} = nothing
)::Nothing
    plot_scalars(
        scalar_name = :strainrate,
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
        plot_type=:CompositionGravity,
        model_output_path=model_output_path,
        material_library_file_path=MATERIAL_COLLECTION.path,
        materials_input_dict=get_materials_input_dict(),
        # General plotting parameters
        dimensions=dimensions,
        xyspacing=xyspacing,
        istart=istart, iend=iend,
        figsize = model_figsize,
        xy_location_contour_legend = (1.0, -100.0), # Contours are not used
        text_box_font_size = 5,
        # General marker parameters
        # (reducing the decimation factor reduces pixelation of the plot)
        marker_size=2.0, decimation_factor=2,
        plot_mesh=0, mesh_line_width=0.1,
        plot_contour_labels=1, contour_line_width=1.0, contour_line_color=:black,
        # Topography parameters
        plot_topography=1, topo_line_width=1.0, 
        topo_line_color=:black,
        # Base level parameters
        plot_base_level=1, base_level_line_width=1.0,
        base_level_line_color=:blue,
        # Plastic strain parameters
        #plot_plastic_strain = true,
        #strain_min = 1.0, strain_max = 6.0,
        #strain_contour_interval = 0.25, strain_cmap=:inferno,
        # Plastic strain rate parameters (log10(1/s))
        plot_plastic_strain_rate = 1,
        strain_rate_min = -18, strain_rate_max = -12,
        strain_rate_contour_interval = 0.5, strain_rate_cmap = :Reds,
        # Sediment age parameters
        plot_sediment_age=1,
        age_min=0.0, age_max=6.0,
        age_contour_interval=0.5, age_cmap=:turbo,
        # Gravity Parameters (mgal)
        gravity_min=-500.0, gravity_max=500.0,
        gravity_spacing=200.0
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
        plot_title="vep_extension_downhill_diffusion",
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
    )
    return nothing
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    Plot.main()
end