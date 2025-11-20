"""
    Plot.jl

Customized plotting module with user defined plotting functions and command-line
plotter that plot output files stored in the model output directory named 
`slab_retreat_output`.

For details on command-line arguments and examples, see the section 
called "Command-line Plotter" or search for `run_cl_plotter` in the EarthBox 
documentation.

# Quick Start

## Plotting via command-line:

```bash
julia Plot.jl <plot_option_name> istart=1 iend=100
```

where `istart` and `iend` are the starting and ending time steps and 
`<plot_option_name>` is the name of the plot option name that can be one of the 
following:

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

include("ReadRun.jl")
import .ReadRun: MATERIAL_MODEL_PATH, MATERIAL_COLLECTION, MODEL_OUTPUT_PATH

const dimensions = (0.0, 1500.0, 0.0, 500.0)
const xyspacing = (100.0, 25.0)
const model_figsize = (10.0, 5.0)

function scalar_plots(;
    model_output_path::String = MODEL_OUTPUT_PATH,
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
end

function velocity_plots(;
    model_output_path::String = MODEL_OUTPUT_PATH,
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
    model_output_path::String = MODEL_OUTPUT_PATH,
    istart::Int64 = 1,
    iend::Union{Int64, Nothing} = nothing
)::Nothing
    plot_markers(
        plot_type=:Composition,
        model_output_path=model_output_path,
        material_library_file_path=MATERIAL_COLLECTION.path,
        material_model_file_path=MATERIAL_MODEL_PATH,
        dimensions=dimensions,
        xyspacing=xyspacing,
        istart=istart, iend=iend,
        # General marker parameters
        marker_size=2.0, decimation_factor=5,
        plot_mesh=0, mesh_line_width=0.1,
        plot_contour_labels=1, contour_line_width=1.0, contour_line_color="black",
        # Topography parameters
        plot_topography=1, topo_line_width=1.0, 
        topo_line_color="red",
        # Temperature contours parameters
        plot_temperature_contours=1,
        temperature_min=100.0, temperature_max=1400.0,
        temperature_contour_interval=100.0,
        temperature_number_format="%6.0f",
        temperature_contour_color="black",
        # Heat flow parameters
        #heatflow_min=0.0, heatflow_max=200.0,
        #heatflow_spacing=50.0,
        #figsize = model_figsize
        # Gravity Parameters (mgal)
        #gravity_min=-1000.0, gravity_max=1000.0,
        #gravity_spacing=200.0
    )
    return nothing
end

function stokes_convergence_plots(
    ;model_output_path::String = MODEL_OUTPUT_PATH,
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
        plot_title="Slab Retreat",
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
        model_output_path = MODEL_OUTPUT_PATH,
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