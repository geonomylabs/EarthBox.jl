module MarkerPlotsManager

include("utils/DataNames.jl")
include("utils/ArrayLookup.jl")
include("utils/PlotMarkerArraysManager.jl")
include("utils/PlotTopoArraysManager.jl")
include("utils/Get.jl")
include("utils/GetModelData.jl")
include("utils/Ticks.jl")
include("utils/JLDData.jl")
include("utils/MarkerColormapManager.jl")
include("filter_plot/FilterPlot.jl")
include("plot_funcs/MarkerPlotFuncs.jl")

import CairoMakie
import EarthBox.ModelDataContainer: ModelData
import EarthBox.UnitConversion: UnitConversionData
import EarthBox.Markers.MarkerMaterials.MaterialsContainer: Materials
import EarthBox.Markers.MarkerMaterials.MaterialsContainer: load_materials_after_checks!
import EarthBox.Markers.MarkerMaterials.MaterialsContainer.MaterialsStateContainer: MaterialsState
import EarthBox.Markers.MarkerMaterials.MaterialsContainer.MaterialsStateContainer: check_state
import EarthBox.ModelDataContainer.Grids2dContainer.ArrayCollection: get_grid_info
import EarthBox.EarthBoxDtypes: MaterialsDictType
import EarthBox.JLDTools: get_jld_data, get_jld_topo_data
import ..PlotParametersManager: PlotParameters
import ..PlotParametersManager: set_parameter_group_attributes!
import ..PlotParametersManager.PlotTimeManager: set_plot_time_info!
import ..PlotParametersManager: set_y_sealevel!
import ..PlotParametersManager: set_base_level_shift!
import ..PlotParametersManager.PlotConversionManager: convert_length_array_units
import ..GridPlotsManager.PlotScalarArraysManager: PlotScalarArrays
import ..GridPlotsManager.PlotScalarArraysManager: set_scalar_arrays!
import ..PlotDtypes: PlotDictType, AxesType
import ..PlotDict: update_plot_dict!
import ..PlotTimeSteppingManager: PlotTimeStepping
import ..PlotTools.InitializePlots: initialize_xy_plot
import ..PlotTools.InitializePlots: initialize_heatflow_composition_plot
import ..PlotTools.InitializePlots: initialize_gravity_composition_plot
import ..PlotTools.InitializePlots: initialize_heatflow_gravity_composition_plot
import ..PlotTools.FinalizePlots: finalize_plot!
import ..PlotTools.PlotUtils: add_text_box
import ..HeatflowPlotsManager
import ..HeatflowPlotsManager: HeatflowPlots
import ..GravityPlotsManager
import ..GravityPlotsManager: GravityPlots
import .PlotMarkerArraysManager: PlotMarkerArrays
import .PlotMarkerArraysManager: set_marker_arrays!
import .MarkerColormapManager: MarkerColorMap
import .PlotTopoArraysManager: PlotTopoArrays
import .DataNames: MarkerDataNames
import .JLDData: get_jld_marker_data, get_jld_topo_filename, get_jld_field_filename
import .GetModelData: get_model_marker_data
import .MarkerPlotFuncs: plot_marker_scalars, plot_temperature_contours
import .MarkerPlotFuncs: plot_topo, plot_base_level, plot_mesh
import .MarkerPlotFuncs: plot_plastic_failure
import .MarkerPlotFuncs: plot_density

Base.@kwdef struct Registry
    Composition::String = "Composition"
    CompositionHeatFlow::String = "CompositionHeatFlow"
    CompositionGravity::String = "CompositionGravity"
    CompositionHeatFlowGravity::String = "CompositionHeatFlowGravity"
    PlasticFailure::String = "PlasticFailure"
    Density::String = "Density"
end

function get_names(registry::Registry)::Vector{String}
    return [
        registry.Composition,
        registry.CompositionHeatFlow,
        registry.CompositionGravity,
        registry.CompositionHeatFlowGravity,
        registry.PlasticFailure,
        registry.Density
    ]
end

mutable struct MarkerPlots
    parameters::PlotParameters
    materials::Materials
    materials_state::MaterialsState
    colormap::MarkerColorMap
    marker_arrays::PlotMarkerArrays
    scalar_arrays::PlotScalarArrays
    topo_arrays::PlotTopoArrays
    unit_conversion::UnitConversionData
    plot_names::Registry
    marker_data_names::MarkerDataNames
    time_stepping::PlotTimeStepping
    heat_flow_plots::HeatflowPlots
    gravity_plots::GravityPlots
    plot_funcs::Dict{String, Function}
end

function MarkerPlots(;
    plot_dict::PlotDictType,
    time_stepping::PlotTimeStepping,
    model_output_path::String,
    plot_output_path::String,
    material_library_file_path::String,
    material_model_file_path::Union{String, Nothing}=nothing,
    materials_input_dict::Union{MaterialsDictType, Nothing}=nothing
)::MarkerPlots
    
    parameters = PlotParameters(plot_dict, model_output_path, plot_output_path)
    materials = Materials()
    materials_state = MaterialsState()
    
    materials_state = load_materials_after_checks!(
        materials, material_library_file_path, material_model_file_path,
        materials_input_dict, materials_state
    )
    
    colormap = MarkerColorMap(materials, materials_state)
    marker_arrays = PlotMarkerArrays()
    scalar_arrays = PlotScalarArrays()
    topo_arrays = PlotTopoArrays()
    unit_conversion = UnitConversionData()
    plot_names = Registry()
    marker_data_names = MarkerDataNames()
    
    heat_flow_plots = HeatflowPlots(
        model_output_path, plot_output_path, time_stepping, marker_data_names)
    
    gravity_plots = GravityPlots(
        model_output_path,
        plot_output_path,
        time_stepping
    )
    
    plot_funcs = make_plot_functions_dict()
    
    return MarkerPlots(
        parameters, materials, materials_state, colormap, marker_arrays,
        scalar_arrays, topo_arrays, unit_conversion, plot_names,
        marker_data_names, time_stepping, heat_flow_plots, gravity_plots,
        plot_funcs
    )
end

function get_keyword_arguments_string()::String
    return """
# Optional Marker Plot Keyword Arguments

## General Marker Plot Parameters
- `marker_size::Float64`: Size of the markers (default: 1.0)
- `plot_contour_labels::Bool`: Plot contour labels (default: false)
- `contour_line_width::Float64`: Width of the contour lines (default: 0.5)
- `contour_line_color::{String, Symbol}`: Contour line color (default: "black")
- `decimation_factor_scatter_overlay::Int`: Decimation factor for scatter overlays (default: 1)
- `nx_contour_grid::Int`: Number of interpolation grid points in x-direction (default: 50)
- `ny_contour_grid::Int`: Number of interpolation grid points in y-direction (default: 50)
- `colorbar_shift_factor::Float64`: Shift factor for color bar (default: 0.0)
- `decimation_factor::Int`: Decimation factor for scatter plots (default: 1)

## Base Level Plot Parameters
- `plot_base_level::Bool`: Whether to plot the base level (default: false)
- `base_level_line_width::Float64`: Width of the base level lines (default: 0.5)
- `base_level_line_color::{String, Symbol}`: Color of the base level lines (default: "blue")

## Topography Plot Parameters
- `plot_topography::Bool`: Whether to plot topography (default: false)
- `topo_line_width::Float64`: Width of the topography lines (default: 0.5)
- `topo_line_color::{String, Symbol}`: Color of the topography lines (default: "black")

## Plastic Strain Plot Parameters
- `plot_plastic_strain::Bool`: Plot marker plastic strain (default: false)
- `strain_min::Float64`: Minimum strain value (default: 1.0)
- `strain_max::Float64`: Maximum strain value (default: 6.0)
- `strain_contour_interval::Float64`: Interval for strain contours (default: 0.25)
- `strain_cmap::{String, Symbol}`: Colormap for strain plot (default: "inferno")

## Plastic Strain Rate Plot Parameters
- `plot_plastic_strain_rate::Bool`: Plot marker plastic strain rate (default: false)
- `strain_rate_min::Float64`: Minimum strain rate value in log10(1/s) (default: -18)
- `strain_rate_max::Float64`: Maximum strain rate value in log10(1/s) (default: -12)
- `strain_rate_contour_interval::Float64`: Interval for strain rate contours in log10(1/s) (default: 1.0)
- `strain_rate_cmap::{String, Symbol}`: Colormap for strain rate plot (default: "Reds")

## Sediment Age Plot Parameters
- `plot_sediment_age::Bool`: Plot marker sediment age (default: false)
- `age_min::Float64`: Minimum age in Myr value (default: 0.0)
- `age_max::Float64`: Maximum age in Myr value (default: 10.0)
- `age_contour_interval::Float64`: Interval for age contours in Myr (default: 1.0)
- `age_cmap::{String, Symbol}`: Colormap for age plot (default: "rainbow")

## Volcanics Age Plot Parameters
- `plot_volcanics_age::Bool`: Plot marker volcanics age (default: false)
- `age_min_volcanics::Float64`: Minimum volcanics age in Myr value (default: 0.0)
- `age_max_volcanics::Float64`: Maximum volcanics age in Myr value (default: 10.0)
- `age_contour_interval_volcanics::Float64`: Interval for volcanics age contours in Myr (default: 1.0)

## Intrusive Age Plot Parameters
- `plot_intrusive_age::Bool`: Plot marker intrusive age (default: false)
- `age_min_intrusive::Float64`: Minimum intrusive age in Myr value (default: 0.0)
- `age_max_intrusive::Float64`: Maximum intrusive age in Myr value (default: 10.0)
- `age_contour_interval_intrusive::Float64`: Interval for intrusive age contours in Myr (default: 1.0)

## Melt Fraction Plot Parameters (only for mantle melting)
- `plot_meltfrac::Bool`: Plot marker melt fraction (default: false)
- `plot_meltfrac_contours::Bool`: Plot melt fraction contours (default: false)

## Extracted Melt Fraction Plot Parameters (only for mantle melting)
- `plot_extracted_meltfrac::Bool`: Plot marker extracted melt fraction (default: false)
- `plot_extracted_meltfrac_contours::Bool`: Plot extracted melt fraction contours (default: false)

## Extractable Melt Fraction Plot Parameters (only for mantle melting)
- `plot_extractable_meltfrac::Bool`: Plot marker extractable melt fraction (default: false)
- `plot_extractable_meltfrac_contours::Bool`: Plot extractable melt fraction contours (default: false)

## Parameters for All Melt-fraction Plots (only for mantle melting)
- `melt_fraction_min::Float64`: Minimum melt fraction value (default: 0.0)
- `melt_fraction_max::Float64`: Maximum melt fraction value (default: 1.0)
- `melt_fraction_contour_interval::Float64`: Interval for melt fraction contours (default: 0.1)
- `meltfrac_cmap::{String, Symbol}`: Colormap for melt fraction plot (default: "inferno")
- `meltfrac_contour_color::{String, Symbol}`: Contour color for melt fraction plot (default: "red")
- `meltfrac_number_format::String`: Number format for melt fraction contours (default: "%6.2f")

## Serpentinization Plot Parameters
- `plot_serpentinization::Bool`: Plot marker serpentinization (default: false)
- `serpentinization_min::Float64`: Minimum serpentinization fraction value (default: 0.0)
- `serpentinization_max::Float64`: Maximum serpentinization fraction value (default: 1.0)
- `serpentinization_contour_interval::Float64`: Interval for serpentinization fraction contours (default: 0.1)
- `serpentinization_cmap::{String, Symbol}`: Colormap for serpentinization plot (default: "inferno")

## Density Plot Parameters
- `plot_density::Bool`: Plot marker density (default: false)
- `density_min::Float64`: Minimum density value in kg/m^3 (default: 0.0)
- `density_max::Float64`: Maximum density value in kg/m^3 (default: 3000.0)
- `density_contour_interval::Float64`: Interval for density contours in kg/m^3 (default: 100.0)

## Grid Temperature Contour Plot Parameters
- `plot_temperature_contours::Bool`: Whether to plot temperature contours (default: false)
- `temperature_min::Float64`: Minimum temperature value in Celsius (default: 0.0)
- `temperature_max::Float64`: Maximum temperature value in Celsius (default: 1300.0)
- `temperature_contour_interval::Float64`: Interval for temperature contours in Celsius (default: 100.0)
- `temperature_contour_color::{String, Symbol}`: Contour color for temperature plot (default: "black")
- `temperature_number_format::String`: Number format for temperature contours (default: "%6.1f")

## Heat flow and gravity subplot parameters
- `heatflow_min::Float64`: Minimum heat flow value in mW/m^2 (default: 0.0)
- `heatflow_max::Float64`: Maximum heat flow value in mW/m^2 (default: 200.0)
- `heatflow_spacing::Float64`: Spacing between heat flow ticks on y-axis in mW/m^2 (default: 50.0)
- `gravity_min::Float64`: Minimum gravity value in mgal (default: -200.0)
- `gravity_max::Float64`: Maximum gravity value in mgal (default: 200.0)
- `gravity_spacing::Float64`: Spacing between gravity ticks on y-axis in mgal (default: 50.0)
- `height_ratios::Vector{Float64}`: Relative height ratios for subplots (default: [0.25, 0.25, 0.75])
"""
end


"""
    plot_markers(
        marker_plots::MarkerPlots;
        plot_type::String="Composition",
        model::Union{ModelData, Nothing}=nothing,
        kwargs...
    )::Nothing

Plot markers for all time Steps.

# Arguments

- `marker_plots::MarkerPlots`:
    - Marker plots object.

- `plot_type::String="Composition"`:
    - Plot type. Must be one of the valid plot types below.

- `model::Union{ModelData, Nothing}=nothing`:
    - Model data object.

$(get_keyword_arguments_string())

# Valid Plot Types
- $(Registry().Composition)
- $(Registry().CompositionHeatFlow)
- $(Registry().CompositionGravity)
- $(Registry().CompositionHeatFlowGravity)
- $(Registry().PlasticFailure)
- $(Registry().Density)

"""
function plot_markers(
    marker_plots::MarkerPlots;
    plot_type::String="Composition",
    model::Union{ModelData, Nothing}=nothing,
    kwargs...
)::Nothing

    println("Plotting markers, plot type: $plot_type")
    
    set_parameters!(marker_plots; kwargs...)
    check_state(marker_plots.materials_state)
    
    # If model is defined then use it for plotting
    if !isnothing(model)
        ioutput = nothing
        call_marker_plot_method!(marker_plots, plot_type, ioutput, model)
    # If model is not defined then use exported model output
    else
        for ioutput in marker_plots.time_stepping.steps
            println(">> Working on marker plot at time step $ioutput")
            call_marker_plot_method!(marker_plots, plot_type, ioutput, model)
        end
    end
    return nothing
end
 
""" 
    set_parameters!(marker_plots::MarkerPlots; kwargs...)::Nothing

Set marker plot parameters.

# Arguments

- marker_plots::MarkerPlots:
    - Marker plots object.

$(get_keyword_arguments_string())

"""
function set_parameters!(marker_plots::MarkerPlots; kwargs...)::Nothing
    update_plot_dict!("marker_plot", marker_plots.parameters.plot_dict; kwargs...)
    set_parameter_group_attributes!(
        marker_plots.parameters.marker_plot_params, 
        marker_plots.parameters.plot_dict["marker_plot"]
        )
end
 
function call_marker_plot_method!(
    marker_plots::MarkerPlots,
    plot_type::String,
    ioutput::Union{Int, Nothing},
    model::Union{ModelData, Nothing}
)::Nothing
    if haskey(marker_plots.plot_funcs, plot_type)
        marker_plots.plot_funcs[plot_type](marker_plots, ioutput, model)
    else
        error("$plot_type is not a valid plot type. Check plot type")
    end
    return nothing
end
 
function make_plot_functions_dict()::Dict{String, Function}
    return Dict(
        "Composition" => make_marker_scalars_plot!,
        "CompositionHeatFlow" => make_marker_scalars_plot_heatflow!,
        "CompositionGravity" => make_marker_scalars_plot_gravity!,
        "CompositionHeatFlowGravity" => make_marker_scalars_plot_heatflow_gravity!,
        "PlasticFailure" => make_plastic_failure_plot!,
        "Density" => make_density_plot!
    )
end
 
 function make_marker_scalars_plot!(
     marker_plots::MarkerPlots,
     ioutput::Union{Int, Nothing},
     model::Union{ModelData, Nothing}
 )::Nothing
    initialize_marker_plot!(marker_plots, ioutput, model)
    fig, axes = initialize_axes_and_mesh!(marker_plots, cmap_name=:blues)
    plot_marker_scalars(
        marker_plots.parameters, marker_plots.marker_arrays,
        marker_plots.materials, marker_plots.colormap, axes
        )
    plot_topo(marker_plots.parameters, marker_plots.topo_arrays, axes)
    plot_base_level(marker_plots.parameters, marker_plots.topo_arrays, axes)
    plot_temperature_contours(marker_plots.parameters, marker_plots.scalar_arrays, axes)
    plot_mesh(marker_plots.parameters, marker_plots.scalar_arrays, axes)
    plot_contour_description!(marker_plots, axes)
    finalize_plot!(
        fig, axes, marker_plots.parameters, "marker_composition",
        extension=marker_plots.parameters.image.extension
        )
    return nothing
 end

 function make_marker_scalars_plot_heatflow!(
     marker_plots::MarkerPlots,
     ioutput::Union{Int, Nothing},
     model::Union{ModelData, Nothing}
 )::Nothing
    initialize_marker_plot!(marker_plots, ioutput, model)
    (
    fig, axes_scatter, axes_heatflow
    ) = initialize_stacked_axes_and_mesh_heatflow!(marker_plots, cmap_name="Blues")
    plot_marker_scalars(
        marker_plots.parameters, marker_plots.marker_arrays, 
        marker_plots.materials, marker_plots.colormap, axes_scatter
    )
    plot_topo(marker_plots.parameters, marker_plots.topo_arrays, axes_scatter)
    plot_base_level(marker_plots.parameters, marker_plots.topo_arrays, axes_scatter)
    plot_temperature_contours(marker_plots.parameters, marker_plots.scalar_arrays, axes_scatter)
    plot_contour_description!(marker_plots, axes_scatter)
    plot_mesh(marker_plots.parameters, marker_plots.scalar_arrays, axes_scatter)
    
    HeatflowPlotsManager.set_ioutput!(marker_plots.heat_flow_plots, ioutput)
    (
    heat_flow_x, hf_gridx, _, _, _, _
    ) = HeatflowPlotsManager.get_heat_flow_grids(marker_plots.heat_flow_plots)
    
    heat_flow_basal_x, _, _, _, _, _ = HeatflowPlotsManager.get_heat_flow_grids(
        marker_plots.heat_flow_plots,
        materials=marker_plots.materials, 
        use_sediment_thickness=true
    )    
    HeatflowPlotsManager.plot_heatflow(
        axes_heatflow, hf_gridx, heat_flow_x, heat_flow_basal_x,
        iplot_basal_heat_flow=false,
        legendfontsize=marker_plots.parameters.fonts.legend_fontsize
    )
     
    finalize_plot!(
        fig, axes_scatter, marker_plots.parameters, "marker_composition_heatflow")
     return nothing
 end
 
function make_marker_scalars_plot_gravity!(
    marker_plots::MarkerPlots,
    ioutput::Union{Int, Nothing},
    model::Union{ModelData, Nothing}
)::Nothing
    initialize_marker_plot!(marker_plots, ioutput, model)
    (
        fig, axes_scatter, axes_gravity
    ) = initialize_stacked_axes_and_mesh_gravity!(marker_plots, cmap_name="Blues")
    plot_marker_scalars(
        marker_plots.parameters, marker_plots.marker_arrays, 
        marker_plots.materials, marker_plots.colormap, axes_scatter
    )
    plot_topo(marker_plots.parameters, marker_plots.topo_arrays, axes_scatter)
    plot_base_level(marker_plots.parameters, marker_plots.topo_arrays, axes_scatter)
    plot_temperature_contours(marker_plots.parameters, marker_plots.scalar_arrays, axes_scatter)
    plot_contour_description!(marker_plots, axes_scatter)
    
    GravityPlotsManager.set_ioutput!(marker_plots.gravity_plots, ioutput)
    (
        gridx, gravity_grid_mgal, gravity_grid_free_air_mgal
    ) = GravityPlotsManager.get_gravity_grids(marker_plots.gravity_plots)

    GravityPlotsManager.plot_gravity(
        axes_gravity, 
        gridx/1000.0,
        gravity_grid_mgal, 
        gravity_grid_free_air_mgal,
        legendfontsize=marker_plots.parameters.fonts.legend_fontsize
    )
    
    finalize_plot!(
        fig, axes_scatter, marker_plots.parameters, "marker_composition_gravity")
    return nothing
end
 
function make_marker_scalars_plot_heatflow_gravity!(
    marker_plots::MarkerPlots,
    ioutput::Union{Int, Nothing},
    model::Union{ModelData, Nothing}
)::Nothing
    initialize_marker_plot!(marker_plots, ioutput, model)
    (
        fig, axes_scatter, axes_heatflow, axes_gravity
    ) = initialize_stacked_axes_and_mesh_heatflow_gravity!(marker_plots, cmap_name="Blues")
    plot_marker_scalars(
        marker_plots.parameters, marker_plots.marker_arrays, 
        marker_plots.materials, marker_plots.colormap, axes_scatter
    )
    plot_topo(marker_plots.parameters, marker_plots.topo_arrays, axes_scatter)
    plot_base_level(marker_plots.parameters, marker_plots.topo_arrays, axes_scatter)
    plot_temperature_contours(marker_plots.parameters, marker_plots.scalar_arrays, axes_scatter)
    plot_contour_description!(marker_plots, axes_scatter)
    
    HeatflowPlotsManager.set_ioutput!(marker_plots.heat_flow_plots, ioutput)
    (
        heat_flow_x, hf_gridx, _, _, _, _
    ) = HeatflowPlotsManager.get_heat_flow_grids(marker_plots.heat_flow_plots)
    heat_flow_basal_x, _, _, _, _, _ = HeatflowPlotsManager.get_heat_flow_grids(
        marker_plots.heat_flow_plots,
        materials=marker_plots.materials, 
        use_sediment_thickness=true
    )
    HeatflowPlotsManager.plot_heatflow(
        axes_heatflow, hf_gridx, heat_flow_x, heat_flow_basal_x,
        legendfontsize=marker_plots.parameters.fonts.legend_fontsize
    )

    GravityPlotsManager.set_ioutput!(marker_plots.gravity_plots, ioutput)
    (
        gridx, gravity_grid_mgal, gravity_grid_free_air_mgal
    ) = GravityPlotsManager.get_gravity_grids(marker_plots.gravity_plots)
    GravityPlotsManager.plot_gravity(
        axes_gravity, gridx/1000.0, gravity_grid_mgal, gravity_grid_free_air_mgal,
        legendfontsize=marker_plots.parameters.fonts.legend_fontsize
    )
    
    finalize_plot!(
        fig, axes_scatter, marker_plots.parameters, "marker_composition_heatflow_gravity")
    return nothing
end
 
 function make_plastic_failure_plot!(
     marker_plots::MarkerPlots,
     ioutput::Union{Int, Nothing},
     model::Union{ModelData, Nothing}
 )::Nothing
    initialize_marker_plot!(marker_plots, ioutput, model)
    axes = initialize_axes_and_mesh!(marker_plots, cmap_name="Blues")
    plot_plastic_failure(
        marker_plots.parameters, marker_plots.marker_arrays, 
        marker_plots.materials, marker_plots.colormap, axes
    )
    plot_temperature_contours(marker_plots.parameters, marker_plots.scalar_arrays, axes)
    finalize_plot!(marker_plots, "marker_plastic", axes)
    return nothing
 end
 
 function make_density_plot!(
     marker_plots::MarkerPlots,
     ioutput::Union{Int, Nothing},
    model::Union{ModelData, Nothing}
 )::Nothing
    initialize_marker_plot!(marker_plots, ioutput, model)
    axes = initialize_axes_and_mesh!(marker_plots, cmap_name="Blues")
    plot_topo!(marker_plots.parameters, marker_plots.topo_arrays, axes)
    plot_density(
        marker_plots.parameters, marker_plots.marker_arrays, 
        marker_plots.materials, axes
    )
    plot_temperature_contours(marker_plots.parameters, marker_plots.scalar_arrays, axes)
    plot_contour_description!(marker_plots, axes)
    finalize_plot!(marker_plots, "marker_density", axes)
    return nothing
 end
 
function initialize_marker_plot!(
    marker_plots::MarkerPlots,
    ioutput::Union{Int, Nothing},
    model::Union{ModelData, Nothing}
)::Nothing
    
    plot_topo = marker_plots.parameters.marker_plot_params.plot_topography
    
    # This option is for loading exported data from the model
    if !isnothing(ioutput)
        update_ioutput!(marker_plots, ioutput)
        load_marker_information_jld!(marker_plots)
        load_grid_information_jld!(marker_plots)
        if plot_topo == 1
            load_topography_information_jld!(marker_plots)
        end
    # This option is for loading data from a ModelData object, which
    # is used for plotting data from a running model
    elseif !isnothing(model)
        ntimestep = model.timestep.parameters.main_time_loop.ntimestep.value
        update_ioutput!(marker_plots, ntimestep)
        load_marker_information_model!(marker_plots, model)
        load_grid_information_model!(marker_plots, model)
        if plot_topo == 1
            load_topography_information_model!(marker_plots, model)
        end
    else
        error("Both ioutput and model cannot be Nothing")
    end
end
 
function initialize_axes_and_mesh!(
    marker_plots::MarkerPlots; 
    cmap_name::Symbol=:blues
)::Tuple{CairoMakie.Figure, CairoMakie.Axis}
    fig, axes = initialize_xy_plot(marker_plots.parameters)
    # Mesh lines must be plotting after scatter to be visible
    plot_mesh(marker_plots.parameters, marker_plots.scalar_arrays, axes)
    return fig, axes
end

function initialize_stacked_axes_and_mesh_heatflow!(
    marker_plots::MarkerPlots; 
    cmap_name::String="Blues"
)::Tuple{CairoMakie.Figure, CairoMakie.Axis, CairoMakie.Axis}
    fig, axes_scatter, axes_heat_flow = initialize_heatflow_composition_plot(marker_plots.parameters)
    # Mesh lines must be plotting after scatter to be visible
    #plot_mesh(marker_plots.parameters, marker_plots.scalar_arrays, axes_scatter)
    return fig, axes_scatter, axes_heat_flow
end

function initialize_stacked_axes_and_mesh_heatflow_gravity!(
    marker_plots::MarkerPlots; 
    cmap_name::String="Blues"
)::Tuple{CairoMakie.Figure, CairoMakie.Axis, CairoMakie.Axis, CairoMakie.Axis}
    (
        fig, axes_scatter, axes_heat_flow, axes_gravity
    ) = initialize_heatflow_gravity_composition_plot(marker_plots.parameters)
    #plot_mesh!(marker_plots.parameters, marker_plots.scalar_arrays, axes_scatter, cmap_name)
    return fig, axes_scatter, axes_heat_flow, axes_gravity
end
 
function initialize_stacked_axes_and_mesh_gravity!(
    marker_plots::MarkerPlots; 
    cmap_name::String="Blues"
)::Tuple{CairoMakie.Figure, CairoMakie.Axis, CairoMakie.Axis}
    fig, axes_scatter, axes_grav = initialize_gravity_composition_plot(marker_plots.parameters)
    #plot_mesh!(marker_plots.parameters, marker_plots.scalar_arrays, axes_scatter, cmap_name)
    return fig, axes_scatter, axes_grav
end

function plot_contour_description!(marker_plots::MarkerPlots, axes::AxesType)::Nothing
    # There is an error in the python module description (ytop and ybottom
    #  are flipped)
    xy_loc = marker_plots.parameters.contours.xy_location_contour_legend
    xloc = xy_loc[1]
    yloc = xy_loc[2]
    description = marker_plots.parameters.contours.contour_description
    # remove the last /n character
    description = description[1:end-1]
    text_box_font_size = marker_plots.parameters.fonts.text_box_font_size
    add_text_box(axes, description, xloc, yloc, fontsize=text_box_font_size)
end

""" Update ioutput plot parameter.

The parameter ioutput is used to access output files for the current
time step if this option is applicable and to define output files.
"""
function update_ioutput!(marker_plots::MarkerPlots, ioutput::Int)::Nothing
    marker_plots.parameters.time.ioutput = ioutput
    return nothing
end
 
""" Load the temperature grid to get mesh information for mesh overlay.

This method also stores temperature data in the scalar_arrays object.
"""
function load_grid_information_jld!(marker_plots::MarkerPlots)::Nothing
    dataname = "TempC"
    jld_filename = get_jld_field_filename(marker_plots.parameters)
    
    _, gridy, gridx, scalar_array, units_dict = get_jld_data(dataname, jld_filename)
    
    set_scalar_arrays!(marker_plots.scalar_arrays, scalar_array)
    
    length_units = units_dict["length_units"]
    
    marker_plots.scalar_arrays.gridx = convert_length_array_units(
        marker_plots.parameters.conversion, length_units, gridx
    )
    
    marker_plots.scalar_arrays.gridy = convert_length_array_units(
        marker_plots.parameters.conversion, length_units, gridy
    )
    return nothing
end
 
function load_grid_information_pymodel!(
    marker_plots::MarkerPlots, 
    model::ModelData
)::Nothing
    
    gridx, gridy, length_units = get_grid_info(model.grids.arrays, "basic")
    
    marker_plots.scalar_arrays.gridx = convert_length_array_units(
        marker_plots.parameters.conversion, length_units, gridx
    )
    
    marker_plots.scalar_arrays.gridy = convert_length_array_units(
        marker_plots.parameters.conversion, length_units, gridy
    )
end

function load_topography_information_jld!(marker_plots::MarkerPlots)::Nothing
    jld_filename = get_jld_topo_filename(marker_plots.parameters)
    topoy, topox, length_units = get_jld_topo_data(jld_filename)
    set_topo_arrays!(marker_plots, topox, topoy, length_units)
end

function load_topography_information_model!(
    marker_plots::MarkerPlots,
    model::ModelData
)::Nothing
    
    topox = copy(model.topography.arrays.gridt.array[1])
    topoy = copy(model.topography.arrays.gridt.array[2])
    length_units = "m"
    set_topo_arrays!(marker_plots, topox, topoy, length_units)
end

function set_topo_arrays!(
    marker_plots::MarkerPlots,
    topox::Vector{Float64},
    topoy::Vector{Float64},
    length_units::String
)::Nothing
    marker_plots.topo_arrays.topox = convert_length_array_units(
        marker_plots.parameters.conversion, length_units, topox
    )
    
    marker_plots.topo_arrays.topoy = convert_length_array_units(
        marker_plots.parameters.conversion, length_units, topoy
    )
    return nothing
end
 
function load_marker_information_jld!(marker_plots::MarkerPlots)::Nothing
    (
        model_time, y_sealevel, base_level_shift, marker_data_dict
    ) = get_jld_marker_data(marker_plots.parameters, marker_plots.marker_data_names)
    set_marker_info!(marker_plots, model_time, y_sealevel, base_level_shift, marker_data_dict)
    return nothing
end
 
function load_marker_information_model!(
    marker_plots::MarkerPlots, 
    model::ModelData
)::Nothing
    (
        model_time, y_sealevel, base_level_shift, marker_data_dict
    ) = get_model_marker_data(
        marker_plots.parameters, marker_plots.marker_data_names, model)
    set_marker_info!(marker_plots, model_time, y_sealevel, base_level_shift, marker_data_dict)
end

function set_marker_info!(
    marker_plots::MarkerPlots,
    model_time::Float64,
    y_sealevel::Float64,
    base_level_shift::Float64,
    marker_data_dict::Dict{String, Vector{Float64}}
)::Nothing
    plot_time_units = marker_plots.parameters.conversion.plot_units.time_units
    set_plot_time_info!(marker_plots.parameters.time, model_time, plot_time_units)
    set_y_sealevel!(marker_plots.parameters, y_sealevel)
    set_base_level_shift!(marker_plots.parameters, base_level_shift)
    set_marker_arrays!(marker_plots.marker_arrays, marker_data_dict)
end

function print_warning_msg(itype_warning::Int, string_info::String)::Nothing
    if itype_warning == 0
        msg = "!!! WARNING !!! $string_info not found in archive."
        println(msg)
    end
end

# """
#     activate_marker_loop_plot!(marker_plots, plot_type="Composition"; kwargs...)
# 
# Activate loop plot for marker plot type.
# """
# function activate_marker_loop_plot!(
#     marker_plots::MarkerPlots,
#     plot_type::String="Composition";
#     kwargs...
# )::Nothing
#     
#     update_plot_dict!("marker_plot", marker_plots.parameters.plot_dict; kwargs...)
#     check_marker_plot_type!(marker_plots, plot_type)
#     set!(marker_plots.parameters.marker_plot_params, 
#          marker_plots.parameters.plot_dict["marker_plot"])
#     check_state!(marker_plots.materials_state)
# end
# 
# """
#     check_marker_plot_type!(marker_plots, plot_type)
# 
# Check if plot_type is a valid marker plot type.
# """
# function check_marker_plot_type!(
#     marker_plots::MarkerPlots, 
#     plot_type::Union{String, Nothing}
# )::Nothing
#     
#     if !isnothing(plot_type)
#         if !(plot_type in get_names(marker_plots.plot_names))
#             error("$plot_type is not a valid marker plot type. Valid options " *
#                   "are $(get_names(marker_plots.plot_names))")
#         end
#     else
#         error("plot_type must be defined for marker plot")
#     end
# end
 
 
end # module
