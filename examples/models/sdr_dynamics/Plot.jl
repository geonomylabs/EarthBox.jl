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
```

Note that `istart` and `iend` are not used by the stokes convergence and yield 
strength plotting functions.

"""
module Plot

using EarthBox

include("Model.jl")
include("Materials.jl")
import .Model: ROOT_PATH_STORAGE
import .Model: xsize, ysize, xo_highres, xf_highres, yf_highres
import .Model: define_case_parameters, PARAMS
import .Materials: get_materials_input_dict, MATERIAL_COLLECTION

#***********************
# User-defined constants
#***********************
#------- Plotting mode constants -------
const use_full_domain = false
const use_gravity_and_heatflow = false # for marker plots only

#------- Miscellaneous constants -------
const xyspacing = use_full_domain ? (25.0, 10.0) : (10.0, 5.0)
const top_padding = 5.0 # km

#------- Plotting output constants -------
const extension = ".png"
const make_additional_pdf_plot = true
const figure_dpi = 200.0
const stflag = use_full_domain ? "zoomout" : "zoomin"

#------- Decimation factors -------
const marker_decimation_factor = use_full_domain ? 5 : 1
const velocity_decimation_factor = use_full_domain ? 8 : 1

#------- Figure sizes -------
# Adjustment factor for vertical size of plots if in full domain mode
const y_factor = use_full_domain ? 1.2 : 1.0
# Take into account the extra space for the gravity and heatflow plots
const model_figsize_markers = use_gravity_and_heatflow ? (12.0, 5.75*y_factor) : (12.0, 3.25*y_factor)
# Scalar plots do not require additional lateral space for color bars
const model_figsize_scalars = (10.0, 3.75*y_factor)

#------- Font sizes -------
const title_fontsize = 20
const axis_title_fontsize = 18
const axis_ticks_fontsize = 16
const colorbar_ticks_fontsize = 16
const colorbar_labels_fontsize = 18
const legend_fontsize = 14
const text_box_font_size = 14

#*****************
# Derived constant
#*****************
const model_center_x = xsize/1000.0/2.0
const plot_width = use_full_domain ? xsize/1000.0 : 100.0 # (xf_highres - xo_highres)/1000.0
const plot_ymin = use_full_domain ? 0.0 : top_padding
const plot_ymax = use_full_domain ? ysize/1000.0 : 40.0 # yf_highres/1000.0

# x-location (km) of the spreading center for each case
const case_spreading_center_x = Dict(
    "case0"  => 200.0,
    "case1"  => 250.0,
    "case2"  => 250.0,
    "case3"  => 250.0,
    "case4"  => 250.0,
    "case5"  => 250.0,
    "case6"  => 250.0,
    "case7"  => 250.0,
    "case8"  => 250.0,
    "case9"  => 250.0,
    "case10" => 250.0,
    "case11" => 250.0,
    "case12" => 250.0,
    "case13" => 250.0,
    "case14" => 250.0,
    "case15" => 250.0,
    "case16" => 250.0,
    "case17" => 250.0,
    "case18" => 250.0,
    "case19" => 250.0,
    "case20" => 250.0,
    "case21" => 250.0,
    "case22" => 250.0,
    "case23" => 250.0,
)

function get_dimensions(case_name::String)::Tuple{Float64, Float64, Float64, Float64}
    spreading_center_x = case_spreading_center_x[case_name]
    xmin = spreading_center_x - plot_width/2.0
    xmax = spreading_center_x + plot_width/2.0
    return (xmin, xmax, plot_ymin, plot_ymax)
end

function scalar_plots(;
    model_output_path::String,
    istart::Int64 = 1,
    iend::Union{Int64, Nothing} = nothing 
)::Nothing
    case_name = get_case_name_from_cl_args()
    plot_dimensions = get_dimensions(case_name)
    plot_scalars(
        scalar_name = :viscosity,
        minimum_value = 18.0,
        maximum_value = 24.0,
        contour_interval = 0.5,
        color_map = "magma",
        stflag = stflag,
        iplot_contour_labels = 0,
        plot_contours = true,
        model_output_path = model_output_path,
        extension=extension, make_pdf=make_additional_pdf_plot,
        dimensions = plot_dimensions,xyspacing = xyspacing,
        istart = istart, iend = iend,
        figsize = model_figsize_scalars, figure_dpi = figure_dpi,
        title_fontsize = title_fontsize, 
        legend_fontsize = legend_fontsize, 
        text_box_font_size = text_box_font_size,
        axis_title_fontsize = axis_title_fontsize, 
        axis_ticks_fontsize = axis_ticks_fontsize,
        colorbar_ticks_fontsize = colorbar_ticks_fontsize,
        colorbar_labels_fontsize = colorbar_labels_fontsize,
    )
    plot_scalars(
        scalar_name = :density,
        minimum_value = 2750.0,
        maximum_value = 3250.0,
        contour_interval = 25.0,
        iplot_contour_labels = 0,
        plot_contours = true,
        color_map = "rainbow",
        stflag = stflag,
        model_output_path = model_output_path,
        extension=extension, make_pdf=make_additional_pdf_plot,
        dimensions = plot_dimensions, xyspacing = xyspacing,
        istart = istart, iend = iend,
        figsize = model_figsize_scalars, figure_dpi = figure_dpi,
        title_fontsize = title_fontsize, 
        legend_fontsize = legend_fontsize, 
        text_box_font_size = text_box_font_size,
        axis_title_fontsize = axis_title_fontsize, 
        axis_ticks_fontsize = axis_ticks_fontsize,
        colorbar_ticks_fontsize = colorbar_ticks_fontsize,
        colorbar_labels_fontsize = colorbar_labels_fontsize,
    )
end

function velocity_plots(;
    model_output_path::String,
    istart::Int64 = 1,
    iend::Union{Int64, Nothing} = nothing
)::Nothing
    plot_velocity(
        model_output_path = model_output_path,
        extension = extension, make_pdf = make_additional_pdf_plot,
        dimensions = plot_dimensions, xyspacing = xyspacing,
        istart = istart, iend = iend,
        decimation_factor = velocity_decimation_factor,
        scale_factor = 10.0,
        figsize = model_figsize_scalars, figure_dpi = figure_dpi,
    )
    return nothing
end

function marker_plots(;
    model_output_path::String,
    istart::Int64 = 1,
    iend::Union{Int64, Nothing} = nothing,
)::Nothing
    # Same case as `run_cl_plotter` / `Model.run_case`: `case_name=` in ARGS (default case0)
    case_name = get_case_name_from_cl_args()
    plot_dimensions = get_dimensions(case_name)

    case_parameters = define_case_parameters(case_name)
    eruption_interval_yr = case_parameters[PARAMS.eruption_interval_yr.name].value
    age_contour_interval_volcanics = eruption_interval_yr/1e6 # Myr
    age_contour_interval_intrusive = age_contour_interval_volcanics

    if use_gravity_and_heatflow
        plot_type = :CompositionHeatFlowGravity
    else
        plot_type = :Composition
    end

    plot_markers(
        plot_type=plot_type,
        model_output_path=model_output_path,
        material_library_file_path=MATERIAL_COLLECTION.path,
        materials_input_dict=get_materials_input_dict(),
        # General plotting parameters
        stflag = stflag,
        extension=extension, make_pdf=make_additional_pdf_plot,
        dimensions=plot_dimensions, xyspacing=xyspacing,
        istart=istart, iend=iend,
        figsize = model_figsize_markers, figure_dpi = figure_dpi,
        xy_location_contour_legend = (-(plot_dimensions[1] + 2.0), -8.0), # Set to negative values to hide
        title_fontsize = 20, axis_title_fontsize = 18, axis_ticks_fontsize = 16,
        colorbar_ticks_fontsize = 16, colorbar_labels_fontsize = 18,
        legend_fontsize = 14, text_box_font_size = 14,
        # General marker parameters
        marker_size=1.0, decimation_factor=marker_decimation_factor,
        plot_mesh=0, mesh_line_width=0.1,
        plot_contour_labels=0, contour_line_width=1.5, contour_line_color="black",
        # Composition colorbar legend customization
        show_composition_matid_labels=false, 
        hidden_composition_matids=[5, 7, 9, 11, 13, 15, 16, 18, 19, 20],
        # Topography parameters
        plot_topography=1, topo_line_width=3.0, 
        topo_line_color="black",
        # Base level parameters
        plot_base_level=0, base_level_line_width=1.0,
        base_level_line_color="blue",
        # Plastic strain parameters
        plot_plastic_strain = true,
        strain_min = 2.0, strain_max = 10.0,
        strain_contour_interval = 0.5, strain_cmap="ice",
        # Plastic strain rate parameters
        plot_plastic_strain_rate = 1,
        strain_rate_min = -18, strain_rate_max = -12,
        strain_rate_contour_interval = 0.5, strain_rate_cmap = "Reds",
        # Intrusive age parameters
        plot_intrusive_age=1,
        age_min_intrusive=0.0, age_max_intrusive=10.0,
        age_contour_interval_intrusive=age_contour_interval_intrusive,
        # Volcanics age parameters
        plot_volcanics_age=1,
        age_min_volcanics=0.0, age_max_volcanics=10.0,
        age_contour_interval_volcanics=age_contour_interval_volcanics,
        # Melt fraction
        plot_meltfrac=1,
        melt_fraction_min=0.0, melt_fraction_max=0.26,
        melt_fraction_contour_interval=0.02,
        use_discontinuous_colormap_meltfrac=true,
        plot_meltfrac_for_gabbro=false,
        meltfrac_cmap="plasma",
        # Temperature contours parameters
        plot_temperature_contours=1,
        temperature_min=100.0, temperature_max=1400.0,
        temperature_contour_interval=100.0,
        temperature_number_format="%6.0f",
        temperature_contour_color="black",
        # Heat flow parameters
        heatflow_min=0.0, heatflow_max=200.0,
        heatflow_spacing=50.0,
        # Gravity Parameters (mgal)
        gravity_min=-400.0, gravity_max=400.0,
        gravity_spacing=100.0
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

function stokes_convergence_plots(
    ;model_output_path::String,
)::Nothing
    plot_stokes_convergence(
        model_output_path=model_output_path,
        figsize=(16, 4), extension=".pdf",
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
        root_path_output = ROOT_PATH_STORAGE,
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
