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

const dimensions = (0.0, 300.0, 0.0, 150.0)
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
        figsize = model_figsize,
        minimum_value = 100.0,
        maximum_value = 1500.0
    )
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
        plot_type=:Composition,
        model_output_path=model_output_path,
        material_library_file_path=MATERIAL_COLLECTION.path,
        materials_input_dict=get_materials_input_dict(),
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
        # Base level parameters
        plot_base_level=1, base_level_line_width=1.0,
        base_level_line_color="blue",
        # Melt fraction parameters
        #plot_meltfrac=1,
        melt_fraction_min=0.0,
        melt_fraction_max=1.0,
        meltfrac_cmap=:turbo, 
        plot_meltfrac_contours=1,
        melt_fraction_contour_interval=0.1,
        meltfrac_contour_color="red", meltfrac_number_format="%6.2f",
        # Temperature contours parameters
        plot_temperature_contours=1,
        temperature_min=100.0, temperature_max=1500.0,
        temperature_contour_interval=100.0,
        temperature_number_format="%6.0f",
        temperature_contour_color="black",
    )
    return nothing
end

function main()::Nothing
    run_cl_plotter(
        root_path_output = ROOT_PATH_OUTPUT,
        marker_plots_func = marker_plots,
        scalar_plots_func = scalar_plots,
        velocity_plots_func = velocity_plots,
    )
    return nothing
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    Plot.main()
end