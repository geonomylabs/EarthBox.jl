module ModelPlots2DManager

include("utils/PlotDtypes.jl")
include("utils/ReadInput.jl")
include("utils/GetIDs.jl")
include("utils/PlotDict.jl")
include("plot_parameters/PlotParametersManager.jl")
include("plot_tools/PlotTools.jl")
include("plot_time_stepping/PlotTimeSteppingManager.jl")
include("grid_plots/GridPlotsManager.jl")
include("heatflow_plots/HeatflowPlotsManager.jl")
include("gravity_plots/GravityPlotsManager.jl")
include("rheology_plots/RheologyPlotsManager.jl")
include("marker_plots/MarkerPlotsManager.jl")
include("plasticity/PlasticityPlotManager.jl")

import EarthBox.PrintFuncs: print_info
import EarthBox.ModelDataContainer: ModelData
import EarthBox.GetPaths: get_model_output_path
import EarthBox.GetArgs: get_case_name_from_cl_args, get_istart, get_iend
import EarthBox.GetArgs: get_model_output_path_from_cl_args, get_plot_option_name
import EarthBox.GetArgs: get_earthbox_project_path_from_args, get_root_path_from_args
import .GridPlotsManager
import .GridPlotsManager: GridPlots
import .GridPlotsManager.ScalarNamesManager: ScalarNames
import .PlotTimeSteppingManager: PlotTimeStepping
import .MarkerPlotsManager
import .MarkerPlotsManager: MarkerPlots
import .RheologyPlotsManager
import .RheologyPlotsManager: RheologyPlots
import .PlasticityPlotManager
import .PlasticityPlotManager: run_plasticity_plotter
import .PlotDict: update_plot_dict!, check_plot_dict, get_plot_dict_template, 
    copy_template_to_plot_dict

function get_general_plotting_keyword_args_string()::String
    return """
# Optional General Plotting Keyword Arguments
- `nsteps::Int`: 
    - Number of time steps to plot
- `nskip::Int:` 
    - Number of time steps to skip
- `istart::Int`: 
    - Initial starting time step
- `iplot_contour_labels::Int`:
    - 0 = off; 1 = plot contour labels
- `color_map::Union{String, Symbol}`:
    - color map name for scalar grid plots (see https://docs.makie.org/dev/explanations/colors)
- `figure_dpi::Float64`:
    - Dots per square inch for figures
- `figsize::Tuple{Float64, Float64}`:
    - Figure size in inches in x-y directions
- `length_units::String`:
    - Length units for plots ("km", "m", "cm")
- `time_units::String`:
    - Time units for plots ("Myr", "yr", "s")
- `velocity_units::String`:
    - Velocity units for plots ("cm/yr", "cm/hr", "m/s", "m/yr", "mm/yr")
- `temperature_units::String`:
    - Temperature units for plots ("C", "K")
- `viscosity_units::String`:
    - Viscosity units for plots ("Pa.s", "log10(Pa.s)")
- `strainrate_units::String`:
    - Strain rate units for plots ("1/s", "log10(1/s)")
- `pressure_units::String`:
    - Pressure units for plots ("Pa", "GPa", "MPa")
- `stress_units::String`:
    - Stress units for plots ("Pa", "GPa", "MPa")
- `linewidth::Float64`:
    - Line width for plots
- `edgecolor::String`:
    - Edge color for plots
- `show_nodes::Bool`:
    - Show nodes on plots
- `dimensions::Tuple{Float64, Float64, Float64, Float64}`:
    - Dimensions of plot: (xmin, xmax, ymin, ymax)
- `dim_close_up::Union{Tuple{Float64, Float64, Float64, Float64}, Nothing}`:
    - Dimensions of zoomed in plot: (xmin, xmax, ymin, ymax)
- `xyspacing::Tuple{Float64, Float64}`:
    - Spacing of ticks for the x and y axes: (xspacing, yspacing)
- `xy_location_contour_legend::Tuple{Float64, Float64}`:
    - Location of low-left-hand corner of contour legend
- `title_fontsize::Int`:
    - Font size for title
- `axis_title_fontsize::Int`:
    - Font size for axis titles
- `axis_labels_fontsize::Int`:
    - Font size for axis labels
- `axis_ticks_fontsize::Int`:
    - Font size for axis ticks
- `contour_label_fontsize::Int`:
    - Font size for contour labels
- `number_format::String`:
    - Number format for contours (e.g. "%6.1f")
- `legend_fontsize::Int`:
    - Font size for legend
- `colorbar_ticks_fontsize::Int`:
    - Font size for color bar ticks
- `colorbar_labels_fontsize::Int`:
    - Font size for color bar labels
- `text_box_font_size::Int`:
    - Font size for text box
- `extension::String`:
    - Extension for plots (e.g. ".png", ".pdf")
"""
end

""" 
    ModelPlots2D

Struct used to manage plot parameters and array data.

# Fields
- `time_stepping::PlotTimeStepping`: 
    - Plot time stepping object
- `grid_plots::GridPlots`: 
    - Grid plots object
- `marker_plots::MarkerPlots`: 
    - Marker plots object
- `rheology_plots::RheologyPlots`: 
    - Rheology plots object
- `scalar_names::ScalarNames`: 
    - Scalar names object
- `scalar_loop_plot_parameters::Dict{String, Dict{String, Union{Float64, Int64, String}}}`: 
    - Dictionary of scalar loop plot parameters
- `active_marker_loop_plot_type::Union{String, Nothing}`: 
    - Active marker loop plot type
- `model_output_path::Union{String, Nothing}`: 
    - Path to model output directory

# Constructor

    ModelPlots2D(;
        plot_output_path::Union{String, Nothing}=nothing,
        material_library_file_path::Union{String, Nothing}=nothing,
        material_model_file_path::Union{String, Nothing}=nothing,
        materials_input_dict::Union{Dict, Nothing}=nothing,
        model_output_path::Union{String, Nothing}=nothing,
        kwargs...
    )

# Arguments
- `plot_output_path::Union{String, Nothing}`:
    - Path of directory where newly created plots will be sent.
- `material_library_file_path::Union{String, Nothing}`:
    - Path of material collection library file with yaml format.
- `material_model_file_path::Union{String, Nothing}`:
    - Path of material input file with yaml format and commonly referred to as materials.yml. 
       If set to nothing then materials_input_dict must be defined.
- `materials_input_dict::Union{Dict, Nothing}`:
    - Dictionary of material inputs where keys are material ID's and each key refers to a 
       dictionary with keys for mat_name, mat_type, mat_domain, red_fraction, green_fraction 
       and blue_fraction
- `model_output_path::Union{String, Nothing}`:
    - Path of directory containing model output files

$(get_general_plotting_keyword_args_string())

"""
struct ModelPlots2D
    time_stepping::PlotTimeStepping
    grid_plots::Union{GridPlots, Nothing}
    marker_plots::Union{MarkerPlots, Nothing}
    rheology_plots::Union{RheologyPlots, Nothing}
    scalar_names::ScalarNames
    scalar_loop_plot_parameters::Dict{String, Dict{String, Union{Float64, Int64, String}}}
    active_marker_loop_plot_type::Union{String, Nothing}
    model_output_path::Union{String, Nothing}
end

function ModelPlots2D(;
    plot_output_path::Union{String, Nothing}=nothing,
    material_library_file_path::Union{String, Nothing}=nothing,
    material_model_file_path::Union{String, Nothing}=nothing,
    materials_input_dict::Union{Dict, Nothing}=nothing,
    model_output_path::Union{String, Nothing}=nothing,
    kwargs...
)
    plot_dict_template = get_plot_dict_template()
    plot_dict = deepcopy(plot_dict_template)
    PlotDict.update_plot_dict!("general_parameters", plot_dict; kwargs...)
    check_plot_dict(plot_dict, plot_dict_template)

    time_stepping = PlotTimeStepping(plot_dict)

    grid_plots = GridPlots(
        plot_dict=plot_dict,
        time_stepping=time_stepping,
        mainpath=model_output_path,
        outpath=plot_output_path
    )

    if !isnothing(material_library_file_path)
        marker_plots = MarkerPlots(
            plot_dict=plot_dict,
            time_stepping=time_stepping,
            model_output_path=model_output_path,
            plot_output_path=plot_output_path,
            material_library_file_path=material_library_file_path,
            material_model_file_path=material_model_file_path,
            materials_input_dict=materials_input_dict
        )
        rheology_plots = RheologyPlots(
            plot_output_path,
            material_library_file_path,
            material_model_file_path=material_model_file_path,
            materials_input_dict=materials_input_dict,
            plot_dict=plot_dict
        )
    else
        marker_plots = nothing
        rheology_plots = nothing
    end

    scalar_names = ScalarNames()
    scalar_loop_plot_parameters = Dict{String, Dict{String, Union{Float64, Int64, String}}}()
    active_marker_loop_plot_type = nothing

    return ModelPlots2D(
        time_stepping,
        grid_plots,
        marker_plots,
        rheology_plots,
        scalar_names,
        scalar_loop_plot_parameters,
        active_marker_loop_plot_type,
        model_output_path
    )
end

function make_stokes_convergence_plot(mp::ModelPlots2D; kwargs...)
    run_plasticity_plotter(
        file_dir=mp.model_output_path;
        kwargs...
    )
end

function get_model_plots_2d(;
    model_output_path::String,
    material_library_file_path::Union{String, Nothing} = nothing,
    material_model_file_path::Union{String, Nothing} = nothing,
    materials_input_dict::Union{Dict, Nothing} = nothing,
    dimensions::Union{Tuple{Float64, Float64, Float64, Float64}, Nothing} = nothing,
    xyspacing::Union{Tuple{Float64, Float64}, Nothing} = nothing,
    istart::Int64 = 1,
    iend::Union{Int64, Nothing} = nothing,
    kwargs...
)::ModelPlots2D
    if isnothing(dimensions)
        dimensions = (0.0, 100.0, 0.0, 100.0)
    end
    if isnothing(xyspacing)
        xyspacing = (10.0, 10.0)
    end
    if isnothing(iend)
        iend = istart
    end
    nsteps = iend - istart + 1
    return ModelPlots2D(
        plot_output_path           = joinpath(model_output_path, "plots"),
        material_library_file_path = material_library_file_path,
        material_model_file_path   = material_model_file_path,
        materials_input_dict       = materials_input_dict,
        model_output_path          = model_output_path,
        nsteps                     = nsteps,
        nskip                      = get(kwargs, :nskip, 0),
        istart                     = istart,
        iplot_contour_labels       = get(kwargs, :iplot_contour_labels, 1),
        color_map                  = get(kwargs, :color_map, "bwr"),
        figure_dpi                 = get(kwargs, :figure_dpi, 150.0),
        figsize                    = get(kwargs, :figsize, (10.0, 4.0)),
        length_units               = get(kwargs, :length_units, "km"),
        time_units                 = get(kwargs, :time_units, "Myr"),
        velocity_units             = get(kwargs, :velocity_units, "cm/yr"),
        temperature_units          = get(kwargs, :temperature_units, "C"),
        viscosity_units            = get(kwargs, :viscosity_units, "log10(Pa.s)"),
        strainrate_units           = get(kwargs, :strainrate_units, "log10(1/s)"),
        pressure_units             = get(kwargs, :pressure_units, "GPa"),
        stress_units               = get(kwargs, :stress_units, "MPa"),
        linewidth                  = get(kwargs, :linewidth, 0.5),
        edgecolor                  = get(kwargs, :edgecolor, "black"),
        show_nodes                 = get(kwargs, :show_nodes, false),
        use_close_up               = get(kwargs, :use_close_up, false),
        dim_close_up               = get(kwargs, :dim_close_up, dimensions),
        dimensions                 = dimensions,
        xyspacing                  = xyspacing,
        xy_location_contour_legend = get(kwargs, :xy_location_contour_legend, (10.0, 10.0)),
        title_fontsize             = get(kwargs, :title_fontsize, 12),
        axis_title_fontsize        = get(kwargs, :axis_title_fontsize, 8),
        axis_labels_fontsize       = get(kwargs, :axis_labels_fontsize, 8),
        axis_ticks_fontsize        = get(kwargs, :axis_ticks_fontsize, 8),
        contour_label_fontsize     = get(kwargs, :contour_label_fontsize, 8),
        number_format              = get(kwargs, :number_format, "%6.1f"),
        legend_fontsize            = get(kwargs, :legend_fontsize, 8),
        colorbar_ticks_fontsize    = get(kwargs, :colorbar_ticks_fontsize, 8),
        colorbar_labels_fontsize   = get(kwargs, :colorbar_labels_fontsize, 8),
        text_box_font_size         = get(kwargs, :text_box_font_size, 8),
        extension                  = get(kwargs, :extension, ".png"),
    )
end

"""
    plot_scalars(;
        scalar_name::Union{String, Symbol},
        model_output_path::String,
        dimensions::Tuple{Float64, Float64, Float64, Float64},
        xy_spacing::Tuple{Float64, Float64},
        istart::Int64 = 1,
        iend::Union{Int64, Nothing} = nothing,
        kwargs...
    )::Nothing

# Arguments
- `scalar_name::Union{String, Symbol}`:
    - Name of scalar to plot
- `model_output_path::String`:
    - Path to output directory containing model output files
- `dimensions::Tuple{Float64, Float64, Float64, Float64}`:
    - Dimensions of plot: (xmin, xmax, ymin, ymax)
- `xyspacing::Tuple{Float64, Float64}`:
    - Spacing of ticks for the x and y axes: (xspacing, yspacing)
- `istart::Int64 = 1`:
    - Initial starting time step
- `iend::Union{Int64, Nothing} = nothing`:
    - Final ending time step

# Optional Scalar Plotting Keyword Arguments
- `contour_interval::Float64`:
    - Interval for contour lines (defaults based on scalar name)
- `minimum_value::Float64`:
    - Minimum value for contour lines (defaults based on scalar name)
- `maximum_value::Float64`:
    - Maximum value for contour lines (defaults based on scalar name)
- `plot_contours::Bool`:
    - Plot contours (default is false)
- `grid_plot_type::String`:
    - Grid plot type: "nomesh" or "mesh" (default is "nomesh")

$(get_general_plotting_keyword_args_string())

# Valid Scalar Names
- `$(GridPlotsManager.ScalarNames().temperature)`
- `$(GridPlotsManager.ScalarNames().viscosity)`
- `$(GridPlotsManager.ScalarNames().strainrate)`
- `$(GridPlotsManager.ScalarNames().pressure)`
- `$(GridPlotsManager.ScalarNames().normal_stress)`
- `$(GridPlotsManager.ScalarNames().shear_stress)`
- `$(GridPlotsManager.ScalarNames().shear_plastic_failure)`
- `$(GridPlotsManager.ScalarNames().normal_plastic_failure)`
- `$(GridPlotsManager.ScalarNames().velocity_x)`
- `$(GridPlotsManager.ScalarNames().velocity_y)`
- `$(GridPlotsManager.ScalarNames().velocity_mag)`
- `$(GridPlotsManager.ScalarNames().density)`
- `$(GridPlotsManager.ScalarNames().thermal_conductivity)`

"""
function plot_scalars(;
    scalar_name::Union{String, Symbol},
    model_output_path::String,
    dimensions::Tuple{Float64, Float64, Float64, Float64},
    xyspacing::Tuple{Float64, Float64},
    istart::Int64 = 1,
    iend::Union{Int64, Nothing} = nothing,
    kwargs...
)::Nothing

    scalar_name = String(scalar_name)
    scalar_names_list = GridPlotsManager.get_scalar_names_list(GridPlotsManager.ScalarNames())
    if !(scalar_name in scalar_names_list)
        error("scalar_name '$(scalar_name)' is not a valid scalar. Valid options are: $(join(scalar_names_list, ", "))")
    end

    println("Plotting scalars for scalar name: $(scalar_name)")

    mp2d = get_model_plots_2d(
        model_output_path          = model_output_path,
        dimensions                 = dimensions, 
        xyspacing                  = xyspacing, 
        istart                     = istart, 
        iend                       = iend;
        kwargs...
    )

    default_scalar_plot_parameters = GridPlotsManager.get_default_scalar_plot_parameters()
    contour_interval_default = default_scalar_plot_parameters[scalar_name].contour_interval
    minimum_value_default = default_scalar_plot_parameters[scalar_name].minimum_value
    maximum_value_default = default_scalar_plot_parameters[scalar_name].maximum_value
    plot_contours_default = default_scalar_plot_parameters[scalar_name].plot_contours
    grid_plot_type_default = default_scalar_plot_parameters[scalar_name].grid_plot_type

    GridPlotsManager.plot_scalars(
        mp2d.grid_plots,
        scalar_name,
        contour_interval = get(kwargs, :contour_interval, contour_interval_default),
        minimum_value    = get(kwargs, :minimum_value, minimum_value_default),
        maximum_value    = get(kwargs, :maximum_value, maximum_value_default),
        plot_contours    = get(kwargs, :plot_contours, plot_contours_default),
        grid_plot_type   = get(kwargs, :grid_plot_type, grid_plot_type_default)
    )
    
end

""" 
    plot_velocity(;
        model_output_path::String,
        dimensions::Tuple{Float64, Float64, Float64, Float64},
        xyspacing::Tuple{Float64, Float64},
        istart::Int64 = 1,
        iend::Union{Int64, Nothing} = nothing,
        kwargs...
    )::Nothing

# Arguments
- `model_output_path::String`:
    - Path to output directory containing model output files
- `dimensions::Tuple{Float64, Float64, Float64, Float64}`:
    - Dimensions of plot: (xmin, xmax, ymin, ymax)
- `xyspacing::Tuple{Float64, Float64}`:
    - Spacing of ticks for the x and y axes: (xspacing, yspacing)
- `istart::Int64 = 1`:
    - Initial starting time step
- `iend::Union{Int64, Nothing} = nothing`:
    - Final ending time step

# Optional Velocity Plotting Keyword Arguments
- `decimation_factor::Int`:
    - Decimation factor for reducing the number of markers plotted
- `scale_factor::Float64`:
    - Scale factor for the velocity plot

$(get_general_plotting_keyword_args_string())

"""
function plot_velocity(;
    model_output_path::String,
    dimensions::Tuple{Float64, Float64, Float64, Float64},
    xyspacing::Tuple{Float64, Float64},
    istart::Int64 = 1,
    iend::Union{Int64, Nothing} = nothing,
    kwargs...
)::Nothing
    mp2d = get_model_plots_2d(
        model_output_path = model_output_path,
        dimensions = dimensions,
        xyspacing = xyspacing,
        istart = istart,
        iend = iend;
        kwargs...
    )
    GridPlotsManager.plot_velocity(
        mp2d.grid_plots,
        decimation_factor=get(kwargs, :decimation_factor, 8),
        scale_factor=get(kwargs, :scale_factor, 10.0)
    )
    return nothing
end

""" 
    plot_markers(;
        plot_type::Union{String, Symbol},
        model_output_path::String,
        material_library_file_path::String,
        material_model_file_path::Union{String, Nothing}=nothing,
        materials_input_dict::Union{Dict, Nothing}=nothing,
        dimensions::Tuple{Float64, Float64, Float64, Float64},
        xyspacing::Tuple{Float64, Float64},
        istart::Int64 = 1,
        iend::Union{Int64, Nothing} = nothing,
        kwargs...
    )::Nothing

# Arguments
- `plot_type::Union{String, Symbol}`:
    - Type of marker plot to plot
- `model_output_path::String`:
    - Path to output directory containing model output files
- `material_library_file_path::String`:
    - Path to material collection library file with yaml format. See 
       [Material Collection Files](@ref) for more information.
- `material_model_file_path::Union{String, Nothing}`:
    - Path to material model input file with yaml format. See 
       [Material Input Files](@ref) for more information.
- `materials_input_dict::Union{Dict, Nothing}`:
    - Dictionary of material inputs where keys are material ID's and each key refers to a 
       dictionary with keys for `mat_name`, `mat_type`, `mat_domain`, 
       `red_fraction`, `green_fraction` and `blue_fraction`. See 
       [Material Input Dictionaries](@ref) for more information.
- `dimensions::Tuple{Float64, Float64, Float64, Float64}`:
    - Dimensions of plot: (xmin, xmax, ymin, ymax)
- `xyspacing::Tuple{Float64, Float64}`:
    - Spacing of ticks for the x and y axes: (xspacing, yspacing)
- `istart::Int64 = 1`:
    - Initial starting time step
- `iend::Union{Int64, Nothing} = nothing`:
    - Final ending time step


# Valid Plot Types
- $(MarkerPlotsManager.Registry().Composition)
    - Marker composition is plotted using RBG values associated with each material
- $(MarkerPlotsManager.Registry().CompositionHeatFlow)
    - Marker composition is plotted using RBG values associated with each material
    - Surface heat flow is plotted above the composition plot.
- $(MarkerPlotsManager.Registry().CompositionGravity)
    - Marker composition is plotted using RBG values associated with each material
    - Surface gravity is plotted above the composition plot.
- $(MarkerPlotsManager.Registry().CompositionHeatFlowGravity)
    - Marker composition is plotted using RBG values associated with each material
    - Surface heat flow and gravity are plotted above the composition plot.
- $(MarkerPlotsManager.Registry().PlasticFailure)
    - Marker plastic failure is plotted using a color map associated with the plastic failure criterion.
- $(MarkerPlotsManager.Registry().Density)
    - Marker density is plotted using a color map associated with the density.


$(MarkerPlotsManager.get_keyword_arguments_string())

$(get_general_plotting_keyword_args_string())

"""
function plot_markers(;
    plot_type::Union{String, Symbol},
    model_output_path::String,
    material_library_file_path::String,
    material_model_file_path::Union{String, Nothing}=nothing,
    materials_input_dict::Union{Dict, Nothing}=nothing,
    dimensions::Tuple{Float64, Float64, Float64, Float64},
    xyspacing::Tuple{Float64, Float64},
    istart::Int64 = 1,
    iend::Union{Int64, Nothing} = nothing,
    kwargs...
)::Nothing

    if !isfile(material_library_file_path)
        error("Material library file does not exist at: $material_library_file_path")
    end
    if isnothing(material_model_file_path) && isnothing(materials_input_dict)
        error("material_model_file_path or materials_input_dict is required")
    end
    if isnothing(materials_input_dict) && !isfile(material_model_file_path)
        error("material_model_file_path does not exist at: $material_model_file_path")
    end

    plot_type = String(plot_type)
    plot_type_names = MarkerPlotsManager.get_names(MarkerPlotsManager.Registry())
    if !(plot_type in plot_type_names)
        error("plot_type '$(plot_type)' is not a valid plot type. Valid options are: $(join(plot_type_names, ", "))")
    end

    mp2d = get_model_plots_2d(
        model_output_path          = model_output_path,
        material_library_file_path = material_library_file_path,
        material_model_file_path   = material_model_file_path,
        materials_input_dict       = materials_input_dict,
        dimensions                 = dimensions,
        xyspacing                  = xyspacing,
        istart                     = istart,
        iend                       = iend;
        kwargs...
    )

    MarkerPlotsManager.plot_markers(
        mp2d.marker_plots, 
        plot_type=plot_type;
        kwargs...
    )
    return nothing

end

"""
    plot_yield_stress(; 
        model_output_path::String,
        kwargs...
    )::Nothing

# Arguments
- `model_output_path::String`:
    - Path to output directory containing model output files

$(RheologyPlotsManager.get_yield_strength_plot_args_string())

"""
function plot_yield_strength(;
    model_output_path::String,
    material_library_file_path::String,
    material_model_file_path::Union{String, Nothing}=nothing,
    materials_input_dict::Union{Dict, Nothing}=nothing,
    kwargs...
)::Nothing
    mp2d = get_model_plots_2d(
        model_output_path = model_output_path,
        material_library_file_path = material_library_file_path,
        material_model_file_path = material_model_file_path,
        materials_input_dict = materials_input_dict; 
        kwargs...
        )
    RheologyPlotsManager.plot_yield_strength(
        mp2d.rheology_plots;
        kwargs...
    )
    return nothing
end

"""
    plot_stokes_convergence(;
        model_output_path::String,
        kwargs...
    )::Nothing

# Arguments
- `model_output_path::String`:
    - Path to output directory containing model output files

$(PlasticityPlotManager.get_stokes_convergence_plot_args_string())

"""
function plot_stokes_convergence(;
    model_output_path::String,
    kwargs...
)::Nothing
    PlasticityPlotManager.run_plasticity_plotter(
        model_output_path;
        kwargs...
    )
    return nothing
end

"""
    run_cl_plotter(;
        model_output_path::Union{String, Nothing}=nothing,
        root_path_output::Union{String, Nothing}=nothing,
        marker_plots_func::Function,
        scalar_plots_func::Function,
        velocity_plots_func::Function,
        stokes_convergence_plots_func::Function,
        yield_strength_plots_func::Function,
    )::Nothing

EarthBox plot runner function that uses command line arguments. Run this function 
from a script that is executed via command line.

# Arguments
- `model_output_path::Union{String, Nothing}`:
    - Path to the model output directory containing the model output files.
       If not provided, then the model output path is defined using the `case_name` 
       and the `root_path` or the model output path provided as a command-line argument.
- `root_path_output::Union{String, Nothing}`:
    - Root path to the model output directory containing the model output files.
       If not provided, then the root path is defined using the command-line 
       argument `root_path=...`.
- `marker_plots_func::Function`:
    - Function to plot markers. 
- `scalar_plots_func::Function`:
    - Function to plot scalars
- `velocity_plots_func::Function`:
    - Function to plot velocity
- `stokes_convergence_plots_func::Function`:
    - Function to plot Stokes convergence
- `yield_strength_plots_func::Function`:
    - Function to plot yield strength
- `model_output_from_script_path::String`:
    - Path to the model output directory sent from the script executing the plotter.
       This is used if the model output path is not specified in the command line arguments.

# Usage 

Assuming the `run_cl_plotter` function is called from a script called Plot.jl ...

To run the plotter with the full model output path specified via a command-line argument:

```bash
julia Plot.jl <plot_option_name> istart=<istart> iend=<iend> model_output_path=<model_output_path> 
```

To run the plotter with the model_output_path defined by `case_name=<case_name>` and
`root_path=<root_path>` from command line arguments:

```bash
julia Plot.jl <plot_option_name> case_name=<case_name> root_path=<root_path> istart=<istart> iend=<iend>
```

To run plotter with `model_output_path` defined by `case_name=<case_name>` from 
command line and the parameter `root_path_output`:

```bash
julia Plot.jl <plot_option_name> case_name=<case_name> istart=<istart> iend=<iend>
```

where:
- `<plot_option_name>`
    - The name of the plot option to use with the following options:
        - `marker_plots`: Plot the markers using the user defined function.
        - `scalar_plots`: Plot the scalars using the user defined function.
        - `velocity_plots`: Plot the velocity using the user defined function.
        - `stokes_convergence_plots`: Plot the Stokes convergence using the user defined function.
        - `yield_strength_plots`: Plot the yield strength using the user defined function.
- `case_name=<case_name>`
    - The name of the model case to plot like `case0`, `case1`, `case2`, etc. If 
       not provided, then the default case name is `case0`.
- `model_output_path=<model_output_path>`
    - The path to the model output directory containing the model output files. 
       If not provided, then the model output path is defined using the `case_name` 
       and the `root_path` or a root path provided as an argument to the 
       `run_cl_plotter` function.
- `root_path=<root_path>`
    - The root path to the model output directory. This will be used along with the 
       `case_name` to define the model output path if the model output path is not 
        provided via command line and the `root_path_output` argument is not 
        defined in the call to the command-line plotter (`run_cl_plotter`).
- `istart=<istart>`
    - The starting time step. Default is 1.
- `iend=<iend>`
    - The ending time step. If not provided, then the ending time step is set to 
       the starting time step.

# Examples

```bash
julia Plot.jl marker_plots istart=1 iend=100 model_output_path=/path/to/naliboff17_case0_output
```

or to use `case_name=<case_name>` and `root_path=<root_path>` to define the 
model output path, run:

```bash
julia Plot.jl marker_plots case_name=case0 root_path=/path/to/root_path istart=1 iend=100
```

or to use `case_name=<case_name>` and the parameter `root_path_output` to define the model 
output path, run:

```bash
julia Plot.jl marker_plots case_name=case0 istart=1 iend=100
```

"""
function run_cl_plotter(;
    model_output_path::Union{String, Nothing}=nothing,
    root_path_output::Union{String, Nothing}=nothing,
    marker_plots_func::Union{Function, Nothing}=nothing,
    scalar_plots_func::Union{Function, Nothing}=nothing,
    velocity_plots_func::Union{Function, Nothing}=nothing,
    stokes_convergence_plots_func::Union{Function, Nothing}=nothing,
    yield_strength_plots_func::Union{Function, Nothing}=nothing,
)::Nothing
    if isnothing(root_path_output)
        root_path_output = get_root_path_from_args()
    end
    if isnothing(model_output_path)
        model_output_path = manage_model_output_path(root_path_output)
    end
    option_name = get_plot_option_name()
    istart = get_istart()
    iend = get_iend()
    print_info("Running command line plotter", level=1)
    print_info("option_name: $option_name", level=2)
    print_info("model_output_path: $model_output_path", level=2)
    print_info("istart: $istart", level=2)
    print_info("iend: $iend", level=2)
    plot_func = get_plot_func(
        option_name,
        marker_plots_func,
        scalar_plots_func,
        velocity_plots_func,
        stokes_convergence_plots_func,
        yield_strength_plots_func
    )
    if option_name == "stokes_convergence_plots" || option_name == "yield_strength_plots"
        if !isnothing(plot_func)
            plot_func(model_output_path = model_output_path)
        else
            println("!!! Warning !!!: plot_func is not defined for option_name: $option_name")
        end
    else
        if !isnothing(plot_func)
            plot_func(model_output_path = model_output_path, istart = istart, iend = iend)
        else
            println("!!! Warning !!!: plot_func is not defined for option_name: $option_name")
        end
    end
    return nothing
end

function manage_model_output_path(root_path_output::String)::String
    # Try to get a case name from the command line arguments. Default to "case0" 
    # if not provided.
    case_name = get_case_name_from_cl_args()
    # Try to get the model output path from the command line arguments. If not
    # provided, then use the case name and the root path output to define the 
    # model output path. Throw errors if the model output path does not exist.
    model_output_path = get_model_output_path(case_name, root_path_output)
    return model_output_path
end

function get_plot_func(
    option_name::String,
    marker_plots_func::Union{Function, Nothing}=nothing,
    scalar_plots_func::Union{Function, Nothing}=nothing,
    velocity_plots_func::Union{Function, Nothing}=nothing,
    stokes_convergence_plots_func::Union{Function, Nothing}=nothing,
    yield_strength_plots_func::Union{Function, Nothing}=nothing
)::Union{Function, Nothing}
    func_dict = Dict(
        "marker_plots" => marker_plots_func,
        "scalar_plots" => scalar_plots_func,
        "velocity_plots" => velocity_plots_func,
        "stokes_convergence_plots" => stokes_convergence_plots_func,
        "yield_strength_plots" => yield_strength_plots_func
    )
    if !haskey(func_dict, option_name)
        throw(ArgumentError(
            "$(option_name) is not a valid option. " *
            "Valid options are $(keys(func_dict))"
        ))
    end
    return func_dict[option_name]
end

end # module