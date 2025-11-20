module GridPlotsManager

include("utils/ScalarNamesManager.jl")
include("utils/DataNames.jl")
include("utils/PlotVelocityArraysManager.jl")
include("utils/PlotScalarArraysManager.jl") # Depends on DataNames
include("utils/PlotVectorArraysManager.jl")
include("utils/PlotGridArraysManager.jl")
include("utils/ArrayLookupManager.jl") # Depends on DataNames
include("velocity_plot/VelocityPlotManager.jl")
include("scalar_plot/ScalarPlotManager.jl")

import EarthBox.ModelDataContainer: ModelData
import ..PlotDtypes: PlotDictType, PlotParametersType
import ..PlotDict: ScalarPlotParameterGroupNames
import .ScalarNamesManager: ScalarNames, get_list as get_scalar_names_list, 
    get_default_scalar_plot_parameters
import .DataNames: ScalarH5Datanames
import .ScalarPlotManager: ScalarPlot
import .ScalarPlotManager: load_params!
import .ScalarPlotManager: make_2dgrid_plots_for_timestep
import .VelocityPlotManager: VelocityPlot
import .VelocityPlotManager: make_2dgrid_plots_for_timestep_velocity
import ..PlotTimeSteppingManager: PlotTimeStepping
import ..PlotParametersManager: update_scalar_plot_parameters!

struct GridPlots
    time_stepping::PlotTimeStepping
    temperature_plot::ScalarPlot
    viscosity_plot::ScalarPlot
    strainrate_plot::ScalarPlot
    pressure_plot::ScalarPlot
    sxx_plot::ScalarPlot
    sxy_plot::ScalarPlot
    plasticflagxy_plot::ScalarPlot
    plasticflagxx_plot::ScalarPlot
    vx_plot::ScalarPlot
    vy_plot::ScalarPlot
    vmag_plot::ScalarPlot
    velocity_plot::VelocityPlot
    density_plot::ScalarPlot
    thermal_conductivity_plot::ScalarPlot
    grid_plot_dict::Dict{String, ScalarPlot}
end

function GridPlots(;
    plot_dict::PlotDictType,
    time_stepping::PlotTimeStepping,
    mainpath::String,
    outpath::String
)::GridPlots
    scalar_names = ScalarNames()
    group_names = ScalarPlotParameterGroupNames()
    h5datanames = ScalarH5Datanames()

    thermal_conductivity_plot = ScalarPlot(
        scalar_name=scalar_names.thermal_conductivity,
        dataname=h5datanames.therm_cond,
        plot_dict=plot_dict,
        parameter_group_name=group_names.thermal_conductivity_plot,
        mainpath=mainpath,
        outpath=outpath
    )

    density_plot = ScalarPlot(
        scalar_name=scalar_names.density,
        dataname=h5datanames.rho_kg_m3,
        plot_dict=plot_dict,
        parameter_group_name=group_names.density_plot,
        mainpath=mainpath,
        outpath=outpath
    )

    temperature_plot = ScalarPlot(
        scalar_name=scalar_names.temperature,
        dataname=h5datanames.TempC,
        plot_dict=plot_dict,
        parameter_group_name=group_names.temperature_plot,
        mainpath=mainpath,
        outpath=outpath
    )

    viscosity_plot = ScalarPlot(
        scalar_name=scalar_names.viscosity,
        dataname=h5datanames.log_eta_Pas,
        plot_dict=plot_dict,
        parameter_group_name=group_names.velocity_plot,
        mainpath=mainpath,
        outpath=outpath
    )

    strainrate_plot = ScalarPlot(
        scalar_name=scalar_names.strainrate,
        dataname=h5datanames.Eii_log,
        plot_dict=plot_dict,
        parameter_group_name=group_names.strain_rate_plot,
        mainpath=mainpath,
        outpath=outpath
    )

    pressure_plot = ScalarPlot(
        scalar_name=scalar_names.pressure,
        dataname=h5datanames.pressure_GPa,
        plot_dict=plot_dict,
        parameter_group_name=group_names.pressure_plot,
        mainpath=mainpath,
        outpath=outpath
    )

    sxx_plot = ScalarPlot(
        scalar_name=scalar_names.normal_stress,
        dataname=h5datanames.Sxx_MPa,
        plot_dict=plot_dict,
        parameter_group_name=group_names.normal_stress_plot,
        mainpath=mainpath,
        outpath=outpath
    )

    sxy_plot = ScalarPlot(
        scalar_name=scalar_names.shear_stress,
        dataname=h5datanames.Sxy_MPa,
        plot_dict=plot_dict,
        parameter_group_name=group_names.shear_stress_plot,
        mainpath=mainpath,
        outpath=outpath
    )

    plasticflagxy_plot = ScalarPlot(
        scalar_name=scalar_names.shear_plastic_failure,
        dataname=h5datanames.plastics,
        plot_dict=plot_dict,
        parameter_group_name=group_names.shear_plastic_failure_plot,
        mainpath=mainpath,
        outpath=outpath
    )

    plasticflagxx_plot = ScalarPlot(
        scalar_name=scalar_names.normal_plastic_failure,
        dataname=h5datanames.plasticn,
        plot_dict=plot_dict,
        parameter_group_name=group_names.normal_plastic_failure_plot,
        mainpath=mainpath,
        outpath=outpath
    )

    vx_plot = ScalarPlot(
        scalar_name=scalar_names.velocity_x,
        dataname=h5datanames.velx_cmyr,
        plot_dict=plot_dict,
        parameter_group_name=group_names.velocity_plot,
        mainpath=mainpath,
        outpath=outpath
    )

    vy_plot = ScalarPlot(
        scalar_name=scalar_names.velocity_y,
        dataname=h5datanames.vely_cmyr,
        plot_dict=plot_dict,
        parameter_group_name=group_names.velocity_plot,
        mainpath=mainpath,
        outpath=outpath
    )

    vmag_plot = ScalarPlot(
        scalar_name=scalar_names.velocity_mag,
        dataname=h5datanames.velmag_cmyr,
        plot_dict=plot_dict,
        parameter_group_name=group_names.velocity_plot,
        mainpath=mainpath,
        outpath=outpath
    )

    velocity_plot = VelocityPlot(
        dataname_x=h5datanames.velx_cmyr,
        dataname_y=h5datanames.vely_cmyr,
        dataname_mag=h5datanames.velmag_cmyr,
        plot_dict=plot_dict,
        parameter_group_name=group_names.velocity_plot,
        mainpath=mainpath,
        outpath=outpath
    )

    grid_plot_dict = Dict(
        temperature_plot.scalar_name => temperature_plot,
        viscosity_plot.scalar_name => viscosity_plot,
        strainrate_plot.scalar_name => strainrate_plot,
        pressure_plot.scalar_name => pressure_plot,
        sxx_plot.scalar_name => sxx_plot,
        sxy_plot.scalar_name => sxy_plot,
        plasticflagxx_plot.scalar_name => plasticflagxx_plot,
        plasticflagxy_plot.scalar_name => plasticflagxy_plot,
        vx_plot.scalar_name => vx_plot,
        vy_plot.scalar_name => vy_plot,
        vmag_plot.scalar_name => vmag_plot,
        density_plot.scalar_name => density_plot,
        thermal_conductivity_plot.scalar_name => thermal_conductivity_plot
    )

    load_all_scalar_plot_parameters!(
        temperature_plot,
        viscosity_plot,
        strainrate_plot,
        pressure_plot,
        sxx_plot,
        sxy_plot,
        plasticflagxy_plot,
        plasticflagxx_plot,
        vx_plot,
        vy_plot,
        vmag_plot,
        density_plot,
        thermal_conductivity_plot,
        plot_dict
    )

    return GridPlots(
        time_stepping,
        temperature_plot,
        viscosity_plot,
        strainrate_plot,
        pressure_plot,
        sxx_plot,
        sxy_plot,
        plasticflagxy_plot,
        plasticflagxx_plot,
        vx_plot,
        vy_plot,
        vmag_plot,
        velocity_plot,
        density_plot,
        thermal_conductivity_plot,
        grid_plot_dict
    )
end

function load_all_scalar_plot_parameters!(
    temperature_plot::ScalarPlot,
    viscosity_plot::ScalarPlot,
    strainrate_plot::ScalarPlot,
    pressure_plot::ScalarPlot,
    sxx_plot::ScalarPlot,
    sxy_plot::ScalarPlot,
    plasticflagxy_plot::ScalarPlot,
    plasticflagxx_plot::ScalarPlot,
    vx_plot::ScalarPlot,
    vy_plot::ScalarPlot,
    vmag_plot::ScalarPlot,
    density_plot::ScalarPlot,
    thermal_conductivity_plot::ScalarPlot,
    plot_dict::PlotDictType
)::Nothing
    load_params!(temperature_plot, plot_dict[temperature_plot.parameter_group_name])
    load_params!(viscosity_plot, plot_dict[viscosity_plot.parameter_group_name])
    load_params!(strainrate_plot, plot_dict[strainrate_plot.parameter_group_name])
    load_params!(pressure_plot, plot_dict[pressure_plot.parameter_group_name])
    load_params!(sxx_plot, plot_dict[sxx_plot.parameter_group_name])
    load_params!(sxy_plot, plot_dict[sxy_plot.parameter_group_name])
    load_params!(plasticflagxy_plot, plot_dict[plasticflagxy_plot.parameter_group_name])
    load_params!(plasticflagxx_plot, plot_dict[plasticflagxx_plot.parameter_group_name])
    load_params!(vx_plot, plot_dict[vx_plot.parameter_group_name])
    load_params!(vy_plot, plot_dict[vy_plot.parameter_group_name])
    load_params!(vmag_plot, plot_dict[vmag_plot.parameter_group_name])
    load_params!(density_plot, plot_dict[density_plot.parameter_group_name])
    load_params!(thermal_conductivity_plot, plot_dict[thermal_conductivity_plot.parameter_group_name])
    return nothing
end

function make_all_gridplots_for_timestep(grid_plots::GridPlots, ioutput::Int)::Nothing
    for (_name, scalar_plot) in grid_plots.grid_plot_dict
        make_2dgrid_plots_for_timestep(scalar_plot; ioutput=ioutput)
    end
    return nothing
end

"""
    plot_scalars(
        grid_plots::GridPlots, 
        scalar_name::String; 
        kwargs...
    )::Nothing

Plot a scalar array for all time Steps.

# Arguments

- `grid_plots::GridPlots`:
    - Grid plots object.
- `scalar_name::String`:
    - Scalar name. Must be one of the valid scalar names below.

# Keyword Arguments

- `model::Union{ModelData, Nothing}=nothing`:
    - Model data object.
- `contour_interval::Float64=1.0`:
    - Contour interval.
- `minimum_value::Float64=0.0`:
    - Minimum value.
- `maximum_value::Float64=100.0`:
    - Maximum value.
- `excluded_value::Float64=-1e38`:
    - Excluded value.
- `plot_contours::Bool=false`:
    - Plot contours.
- `grid_plot_type::String="nomesh"`:
    - Grid plot type: "nomesh" or "mesh"

# Valid Scalar Names
- `$(ScalarNames().temperature)`
- `$(ScalarNames().viscosity)`
- `$(ScalarNames().strainrate)`
- `$(ScalarNames().pressure)`
- `$(ScalarNames().normal_stress)`
- `$(ScalarNames().shear_stress)`
- `$(ScalarNames().shear_plastic_failure)`
- `$(ScalarNames().normal_plastic_failure)`
- `$(ScalarNames().velocity_x)`
- `$(ScalarNames().velocity_y)`
- `$(ScalarNames().velocity_mag)`
- `$(ScalarNames().density)`
- `$(ScalarNames().thermal_conductivity)`

"""
function plot_scalars(
    grid_plots::GridPlots,
    scalar_name::String;
    model::Union{ModelData, Nothing}=nothing,
    contour_interval::Float64=1.0,
    minimum_value::Float64=0.0,
    maximum_value::Float64=100.0,
    excluded_value::Float64=-1e38,
    plot_contours::Bool=false,
    grid_plot_type::String="nomesh"
)::Nothing
    plot_parameters = get_plot_parameters_dict(
        contour_interval,
        minimum_value,
        maximum_value,
        excluded_value,
        plot_contours,
        grid_plot_type
    )

    if isnothing(model)
        plot_output_scalars(grid_plots, scalar_name, plot_parameters)
    else
        make_scalar_plot(
            grid_plots,
            scalar_name,
            plot_parameters,
            model=model
        )
    end
    return nothing
end

function get_plot_parameters_dict(
    contour_interval::Float64,
    minimum_value::Float64,
    maximum_value::Float64,
    excluded_value::Float64,
    plot_contours::Bool,
    grid_plot_type::String
)::PlotParametersType
    return Dict(
        "contour_interval" => contour_interval,
        "minimum_value" => minimum_value,
        "maximum_value" => maximum_value,
        "excluded_value" => excluded_value,
        "iplot_contours" => get_integer_flag(plot_contours),
        "grid_plot_type" => grid_plot_type
    )
end

function plot_output_scalars(
    grid_plots::GridPlots,
    scalar_name::String,
    plot_parameters::PlotParametersType
)::Nothing
    for ioutput in grid_plots.time_stepping.steps
        make_scalar_plot(grid_plots, scalar_name, plot_parameters, ioutput=ioutput)
    end
    return nothing
end

function get_integer_flag(option::Union{Bool, Nothing})::Union{Int64, Nothing}
    if !isnothing(option)
        iflag = 0
        if option
            iflag = 1
        end
        return iflag
    end
    return nothing
end

function make_scalar_plot(
    grid_plots::GridPlots,
    scalar_name::String,
    plot_parameters::PlotParametersType;
    ioutput::Union{Int64, Nothing}=nothing,
    model::Union{ModelData, Nothing}=nothing
)::Nothing
    grid_plot_name = get_grid_plot_name(grid_plots, scalar_name)
    grid_plot = grid_plots.grid_plot_dict[grid_plot_name]
    grid_plot.parameters.options.iplot = 1
    update_scalar_plot_parameters!(grid_plot.parameters, plot_parameters)
    make_2dgrid_plots_for_timestep(grid_plot, ioutput=ioutput, model=model)
    return nothing
end

function get_grid_plot_name(grid_plots::GridPlots, scalar_name::String)::String
    scalar_names = ScalarNames()
    scalar_name_to_grid_plot_name = Dict(
        scalar_names.temperature => grid_plots.temperature_plot.scalar_name,
        scalar_names.viscosity => grid_plots.viscosity_plot.scalar_name,
        scalar_names.strainrate => grid_plots.strainrate_plot.scalar_name,
        scalar_names.pressure => grid_plots.pressure_plot.scalar_name,
        scalar_names.normal_stress => grid_plots.sxx_plot.scalar_name,
        scalar_names.shear_stress => grid_plots.sxy_plot.scalar_name,
        scalar_names.shear_plastic_failure => grid_plots.plasticflagxy_plot.scalar_name,
        scalar_names.normal_plastic_failure => grid_plots.plasticflagxx_plot.scalar_name,
        scalar_names.velocity_x => grid_plots.vx_plot.scalar_name,
        scalar_names.velocity_y => grid_plots.vy_plot.scalar_name,
        scalar_names.velocity_mag => grid_plots.vmag_plot.scalar_name,
        scalar_names.density => grid_plots.density_plot.scalar_name,
        scalar_names.thermal_conductivity => grid_plots.thermal_conductivity_plot.scalar_name
    )
    return scalar_name_to_grid_plot_name[scalar_name]
end

"""
    plot_velocity(grid_plots::GridPlots; kwargs...)::Nothing

# Arguments

- `grid_plots::GridPlots`:
    - Grid plots object.

# Keyword Arguments

- `decimation_factor::Int64` = 1:
    - Decimation factor for velocity vectors.
- `scale_factor::Float64` = 10.0:
    - Scale factor for velocity vectors.

Plot the velocity vectors.
"""
function plot_velocity(
    grid_plots::GridPlots;
    decimation_factor::Int64=1,
    scale_factor::Float64=10.0
)::Nothing
    println("Working on plots for velocity")
    for ioutput in grid_plots.time_stepping.steps
        make_vector_plot(
            grid_plots,
            ioutput,
            decimation_factor=decimation_factor,
            scale_factor=scale_factor
        )
    end
    return nothing
end

function make_vector_plot(
    grid_plots::GridPlots,
    ioutput::Int64;
    decimation_factor::Int64=1,
    scale_factor::Float64=10.0
)::Nothing
    make_2dgrid_plots_for_timestep_velocity(
        velocity_plot=grid_plots.velocity_plot,
        ioutput=ioutput,
        decimation_factor=decimation_factor,
        scale_factor=scale_factor
    )
    return nothing
end

end # module 