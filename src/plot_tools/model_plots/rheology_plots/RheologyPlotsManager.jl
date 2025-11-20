module RheologyPlotsManager

include("utils/GetBasicLithosphereIDs.jl")
include("utils/GetRheology.jl")
include("utils/CompositeYieldStress.jl")

import CairoMakie
import EarthBox.Markers.MarkerMaterials.MaterialsContainer: Materials
import EarthBox.Markers.MarkerMaterials.MaterialsContainer: load_materials_after_checks!
import EarthBox.Markers.MarkerMaterials.MaterialsContainer.MaterialsStateContainer: MaterialsState
import EarthBox.Markers.MarkerMaterials.MaterialsContainer.MaterialsStateContainer: check_state
import EarthBox.EarthBoxDtypes: MaterialsDictType
import EarthBox.SurfaceProcesses.Sealevel.RelativeBaseLevel.ReferenceLithosphere: 
    make_lithosphere_model_from_user_inputs
import ..PlotParametersManager: PlotParameters
import ..PlotDtypes: PlotDictType, AxesType
import ..PlotTools.InitializePlots: initialize_xy_plot
import ..PlotTools.FinalizePlots: finalize_plot!
import .CompositeYieldStress: calculate_lithosphere_composite_yield_stress
import .GetBasicLithosphereIDs: get_lithosphere_id_tuple, get_material_id_grid
import .GetRheology: get_basic_lithosphere_rheology

""" Rheology plots struct for managing rheology plots.

# Fields
- `parameters::PlotParameters`: Plot parameters object
- `materials::Materials`: Materials object
- `materials_state::MaterialsState`: Materials state object
- `plot_output_path::String`: Path to output directory for plots
"""
mutable struct RheologyPlots
    parameters::PlotParameters
    materials::Materials
    materials_state::MaterialsState
    plot_output_path::String
end

""" Constructor for RheologyPlots.

# Arguments
- `plot_output_path::String`: Path to output directory for plots
- `material_library_file_path::String`: Path to material library file
- `material_model_file_path::Union{String, Nothing}`: Path to material model file
- `materials_input_dict::Union{MaterialsDictType, Nothing}`: Materials input dictionary
- `plot_dict::PlotDictType`: Plot dictionary with parameters
"""
function RheologyPlots(
    plot_output_path::String,
    material_library_file_path::String;
    material_model_file_path::Union{String, Nothing}=nothing,
    materials_input_dict::Union{MaterialsDictType, Nothing}=nothing,
    plot_dict::PlotDictType=Dict{String, Dict{String, Union{Float64, Int64, String, Bool, Tuple{Float64, Vararg{Float64}}}}}()
)::RheologyPlots
    
    parameters = PlotParameters(plot_dict, "", plot_output_path)
    materials = Materials()
    materials_state = MaterialsState()
    
    materials_state = load_materials_after_checks!(
        materials, material_library_file_path, material_model_file_path,
        materials_input_dict, materials_state
    )
    
    return RheologyPlots(parameters, materials, materials_state, plot_output_path)
end

function get_yield_strength_plot_args_string()::String
    return """
# Yield Strength Plot Keyword Arguments

## Plot Parameters
- `figsize::Tuple{Float64, Float64}`: Figure size (default: (5.0, 10.0))
- `depth_plot_limit::Float64`: Depth plot limit in km (default: 160.0)
- `depth_axis_spacing::Float64`: Depth axis spacing in km (default: 10.0)
- `temperature_plot_limit::Float64`: Temperature plot limit in Celsius (default: 1400.0)
- `temperature_axis_spacing::Float64`: Temperature axis spacing in Celsius (default: 100.0)
- `yield_stress_plot_limit::Float64`: Yield stress plot limit in MPa (default: 500.0)
- `yield_stress_axis_spacing::Float64`: Yield stress axis spacing in MPa (default: 50.0)

## Strain Rate Parameters
- `strain_rate::Float64`: Strain rate invariant in 1/s (default: 1.0e-15)

## Geotherm Model Option
- `iuse_linear_segments::Bool`: Use a temperature profile with four linear segments. 
  If false, an analytical three-layer temperature profile is used instead.

## Thickness Parameters
- `thickness_upr_cont_crust_meters::Float64`: Thickness of upper continental crust in meters
- `thickness_lwr_cont_crust_meters::Float64`: Thickness of lower continental crust in meters
- `thickness_lithosphere_meters::Float64`: Thickness of lithosphere in meters
- `thickness_thermal_lithosphere::Float64`: Thickness of thermal lithosphere in meters
- `thickness_asthenosphere_meters::Float64`: Thickness of asthenosphere in meters
- `dy_meters::Float64`: Grid spacing in meters

## Density Parameters
- `expansivity::Float64`: Thermal expansivity in 1/K
- `compressibility::Float64`: Compressibility in 1/Pa
- `density_upper_continental_crust::Float64`: Density of upper continental crust in kg/m^3
- `density_lower_continental_crust::Float64`: Density of lower continental crust in kg/m^3
- `density_mantle_lithosphere::Float64`: Density of mantle lithosphere in kg/m^3
- `density_asthenosphere::Float64`: Density of asthenosphere in kg/m^3

## Temperature Parameters
- `temperature_top_celsius::Float64`: Temperature at top of lithosphere in Celsius
- `temperature_moho_celsius::Float64`: Temperature at Moho in Celsius
- `temperature_base_lith_celsius::Float64`: Temperature at base of lithosphere in Celsius
- `adiabatic_gradient_kelvin_km::Float64`: Adiabatic gradient in Kelvin/km

## Thermal Conductivity (Analytical 3-layer Model)
- `conductivity_upper_crust::Float64`: Thermal conductivity of upper crust in W/m/K
- `conductivity_lower_crust::Float64`: Thermal conductivity of lower crust in W/m/K
- `conductivity_mantle::Float64`: Thermal conductivity of mantle in W/m/K

## Heat Production (Analytical 3-layer Model)
- `heat_production_upper_crust::Float64`: Heat production in upper crust in W/m^3
- `heat_production_lower_crust::Float64`: Heat production in lower crust in W/m^3
- `heat_production_mantle::Float64`: Heat production in mantle in W/m^3

## Fluid Pressure Option
- `iuse_fluid_pressure_for_yield::Int`: Use fluid pressure in yield stress calculation (default: 0)
"""
end

""" 
    plot_yield_strength(rheology_plots::RheologyPlots; kwargs...)::Nothing

Plot yield strength profile.

# Arguments

- `rheology_plots::RheologyPlots`:
    - Rheology plots object.

$(get_yield_strength_plot_args_string())

"""
function plot_yield_strength(
    rheology_plots::RheologyPlots;
    kwargs...
)::Nothing
    
    check_kwargs(; kwargs...)
    
    lith_ids = get_lithosphere_id_tuple(rheology_plots.materials)
    
    (
        upr_crust_rheology,
        lwr_crust_rheology,
        mantle_lithosphere_rheology,
        asthenosphere_rheology
    ) = get_basic_lithosphere_rheology(rheology_plots.materials, lith_ids)
    
    (
        gridy, temp_gridy,
        _density_gridy, pressure_gridy,
        lithosphere_thicknesses
    ) = make_lithosphere_model_from_user_inputs(; kwargs...)
    
    conductivity_upper_crust = get(kwargs, :conductivity_upper_crust, 2.5)
    # Top temperature gradient
    dy = gridy[2] - gridy[1]
    temp_gradient = (temp_gridy[2] - temp_gridy[1]) / dy
    surface_heat_flow = conductivity_upper_crust * temp_gradient * 1000.0
    
    matid_gridy = get_material_id_grid(gridy, lithosphere_thicknesses, lith_ids)
    
    strain_rate = get(kwargs, :strain_rate, 1.0e-15)

    figure_dpi = rheology_plots.parameters.image.figure_dpi
    figsize = get(kwargs, :figsize, (5.0, 10.0))
    figsize_pixels = (figsize[1] * figure_dpi, figsize[2] * figure_dpi)
    
    (
        yield_stress_gridy_no_damage,
        dislocation_stress,
        diffusion_stress,
        peierls_stress
    ) = calculate_lithosphere_composite_yield_stress(
        strain_rate,
        lith_ids,
        matid_gridy,
        upr_crust_rheology,
        lwr_crust_rheology,
        mantle_lithosphere_rheology,
        asthenosphere_rheology,
        temp_gridy,
        pressure_gridy,
        gridy;
        use_damaged_state=false,
        iuse_fluid_pressure_for_yield=get(kwargs, :iuse_fluid_pressure_for_yield, 0)
    )
    
    (
        yield_stress_gridy_with_damage,
        _, _, _
    ) = calculate_lithosphere_composite_yield_stress(
        strain_rate,
        lith_ids,
        matid_gridy,
        upr_crust_rheology,
        lwr_crust_rheology,
        mantle_lithosphere_rheology,
        asthenosphere_rheology,
        temp_gridy,
        pressure_gridy,
        gridy;
        use_damaged_state=true,
        iuse_fluid_pressure_for_yield=get(kwargs, :iuse_fluid_pressure_for_yield, 0)
    )
    
    plastic_failure = upr_crust_rheology.plastic_failure
    cohesions = [
        plastic_failure.cohesion_initial,
        plastic_failure.cohesion_final
    ]
    
    friction_angles = [
        plastic_failure.friction_angle_initial,
        plastic_failure.friction_angle_final
    ]
    
    check_direc(rheology_plots.plot_output_path)
    make_strength_profile_plot_twin(
        rheology_plots.plot_output_path,
        gridy,
        yield_stress_gridy_no_damage,
        yield_stress_gridy_with_damage,
        dislocation_stress,
        diffusion_stress,
        peierls_stress,
        temp_gridy,
        cohesions,
        friction_angles,
        surface_heat_flow,
        figsize_pixels;
        kwargs...
    )
    
    make_effective_viscosity_profile_plot(
        rheology_plots.plot_output_path,
        gridy,
        yield_stress_gridy_no_damage,
        yield_stress_gridy_with_damage,
        cohesions,
        friction_angles,
        figsize_pixels;
        kwargs...
    )
    
    return nothing
end

""" Check kwargs for valid keys.

# Arguments
- `kwargs...`: Keyword arguments to check
"""
function check_kwargs(; kwargs...)::Nothing
    valid_keys = [
        :figsize,
        :depth_plot_limit,
        :depth_axis_spacing,
        :temperature_plot_limit,
        :temperature_axis_spacing,
        :yield_stress_plot_limit,
        :yield_stress_axis_spacing,
        :log10_viscosity_plot_min,
        :log10_viscosity_plot_max,
        :log10_viscosity_axis_spacing,
        :strain_rate,
        :thickness_upr_cont_crust_meters,
        :thickness_lwr_cont_crust_meters,
        :thickness_lithosphere_meters,
        :thickness_asthenosphere_meters,
        :dy_meters,
        :expansivity,
        :compressibility,
        :density_upper_continental_crust,
        :density_lower_continental_crust,
        :density_mantle_lithosphere,
        :density_asthenosphere,
        :temperature_top_celsius,
        :temperature_moho_celsius,
        :temperature_base_lith_celsius,
        :adiabatic_gradient_kelvin_km,
        :conductivity_upper_crust,
        :conductivity_lower_crust,
        :conductivity_mantle,
        :heat_production_upper_crust,
        :heat_production_lower_crust,
        :heat_production_mantle,
        :thickness_thermal_lithosphere,
        :iuse_linear_segments,
        :plot_deformation_mechanisms,
        :iuse_fluid_pressure_for_yield,
        :extension
    ]
    
    for (key, _) in kwargs
        if !(key in valid_keys)
            println("Valid inputs for kwargs are:")
            println(valid_keys)
            error("Invalid key: $key. See the list of valid keys above.")
        end
    end
end

""" Make strength profile plot with different x-axes.

# Arguments
- `plot_output_path::String`: Path to output directory
- `gridy::Vector{Float64}`: Y-coordinate grid
- `yield_stress_gridy_no_damage::Vector{Float64}`: Yield stress without damage
- `yield_stress_gridy_with_damage::Vector{Float64}`: Yield stress with damage
- `dislocation_stress::Vector{Float64}`: Dislocation stress
- `diffusion_stress::Vector{Float64}`: Diffusion stress
- `_peierls_stress::Vector{Float64}`: Peierls stress (unused)
- `temp_gridy::Vector{Float64}`: Temperature grid
- `cohesions::Vector{Float64}`: Cohesion values
- `friction_angles::Vector{Float64}`: Friction angle values
- `surface_heat_flow::Float64`: Surface heat flow
- `kwargs...`: Additional keyword arguments
"""
function make_strength_profile_plot_twin(
    plot_output_path::String,
    gridy::Vector{Float64},
    yield_stress_gridy_no_damage::Vector{Float64},
    yield_stress_gridy_with_damage::Vector{Float64},
    dislocation_stress::Vector{Float64},
    diffusion_stress::Vector{Float64},
    _peierls_stress::Vector{Float64},
    temp_gridy::Vector{Float64},
    cohesions::Vector{Float64},
    friction_angles::Vector{Float64},
    surface_heat_flow::Float64,
    figsize_pixels::Tuple{Float64, Float64};
    kwargs...
)::Nothing
    
    depth_plot_limit = get(kwargs, :depth_plot_limit, 160.0)
    depth_axis_spacing = get(kwargs, :depth_axis_spacing, 10.0)
    temperature_plot_limit = get(kwargs, :temperature_plot_limit, 1400.0)
    temperature_axis_spacing = get(kwargs, :temperature_axis_spacing, 100.0)
    yield_stress_plot_limit = get(kwargs, :yield_stress_plot_limit, 500.0)
    yield_stress_axis_spacing = get(kwargs, :yield_stress_axis_spacing, 50.0)
    plot_deformation_mechanisms = get(kwargs, :plot_deformation_mechanisms, false)
    extension = get(kwargs, :extension, ".png")
    #figsize = get(kwargs, :figsize, (5.0, 10.0))
    
    fig, ax1 = initialize_plot(figsize_pixels)
    
    # Convert Pa to MPa
    yield_stress_no_damage_mpa = yield_stress_gridy_no_damage ./ 1e6
    yield_stress_with_damage_mpa = yield_stress_gridy_with_damage ./ 1e6
    dislocation_stress_mpa = dislocation_stress ./ 1e6
    diffusion_stress_mpa = diffusion_stress ./ 1e6
    
    # Convert m to km
    depth_km = gridy ./ 1000.0
    
    # Plot yield stress profiles
    CairoMakie.lines!(ax1, yield_stress_no_damage_mpa, depth_km, 
           label="Yield Stress (No Damage)", color=:black)
    CairoMakie.lines!(ax1, yield_stress_with_damage_mpa, depth_km, 
           label="Yield Stress (Damaged)", color=:black, linestyle=:dot)
    
    if plot_deformation_mechanisms
        CairoMakie.lines!(ax1, dislocation_stress_mpa, depth_km, 
               label="Dislocation Stress", color=:orange, linestyle=:dash)
        CairoMakie.lines!(ax1, diffusion_stress_mpa, depth_km, 
               label="Diffusion Stress", color=:green, linestyle=:dashdot)
    end
    
    ax1.xlabel = "Yield Stress (MPa)"
    ax1.ylabel = "Depth (km)"
    CairoMakie.ylims!(ax1, 0, depth_plot_limit)
    ax1.yticks = get_ticks(depth_plot_limit, depth_axis_spacing)
    ax1.yreversed = true
    CairoMakie.xlims!(ax1, 0, yield_stress_plot_limit)
    ax1.xticks = get_ticks(yield_stress_plot_limit, yield_stress_axis_spacing)
    CairoMakie.axislegend(ax1, position=:rt, framevisible=false, labelsize=10)
    
    ax2 = CairoMakie.Axis(fig[1, 1], xaxisposition=:top, yaxisposition=:right)
    CairoMakie.hideydecorations!(ax2)
    temp_celsius = temp_gridy .- 273.15
    CairoMakie.lines!(ax2, temp_celsius, depth_km, label="Temperature", 
           color=:red, linestyle=:dash, )
    ax2.xlabel = "Temperature (Celsius)"
    CairoMakie.xlims!(ax2, 0, temperature_plot_limit)
    ax2.xticks = get_ticks(temperature_plot_limit, temperature_axis_spacing)
    CairoMakie.ylims!(ax2, 0, depth_plot_limit)
    ax2.yticks = get_ticks(depth_plot_limit, depth_axis_spacing)
    ax2.yreversed = true
    ax2.xgridvisible = true
    ax2.xgridstyle = :dash
    ax2.xminorgridvisible = true
    ax2.xminorticks = CairoMakie.IntervalsBetween(2)
    CairoMakie.axislegend(ax2, position=(0.78, 0.94), framevisible=false, labelsize=10)
    
    strain_rate = get(kwargs, :strain_rate, 1.0e-15)
    text_str = (
        "Strain Rate: $strain_rate 1/s\n" *
        "Cohesion (No Damage): $(cohesions[1]/1e6) MPa\n" *
        "Cohesion (Damaged): $(cohesions[2]/1e6) MPa\n" *
        "Friction Angle (No Damage): $(friction_angles[1]) degrees\n" *
        "Friction Angle (Damaged): $(friction_angles[2]) degrees\n" *
        "Surface Heat Flow: $(round(surface_heat_flow, digits=2)) mW/m^2"
    )
    
    CairoMakie.text!(ax1, text_str, position=(0.02, 0.02), 
          fontsize=8, align=(:left, :bottom))
    
    plot_name = joinpath(plot_output_path, "yield_stress" * extension)
    CairoMakie.save(plot_name, fig)
    #CairoMakie.close(fig)
    
    return nothing
end

""" Make effective viscosity profile plot.

# Arguments
- `plot_output_path::String`: Path to output directory
- `gridy::Vector{Float64}`: Y-coordinate grid
- `yield_stress_gridy_no_damage::Vector{Float64}`: Yield stress without damage
- `yield_stress_gridy_with_damage::Vector{Float64}`: Yield stress with damage
- `cohesions::Vector{Float64}`: Cohesion values
- `friction_angles::Vector{Float64}`: Friction angle values
- `kwargs...`: Additional keyword arguments
"""
function make_effective_viscosity_profile_plot(
    plot_output_path::String,
    gridy::Vector{Float64},
    yield_stress_gridy_no_damage::Vector{Float64},
    yield_stress_gridy_with_damage::Vector{Float64},
    cohesions::Vector{Float64},
    friction_angles::Vector{Float64},
    figsize_pixels::Tuple{Float64, Float64};
    kwargs...
)::Nothing
    
    strain_rate = get(kwargs, :strain_rate, 1.0e-15)
    depth_plot_limit = get(kwargs, :depth_plot_limit, 160.0)
    depth_axis_spacing = get(kwargs, :depth_axis_spacing, 10.0)
    log10_viscosity_plot_min = get(kwargs, :log10_viscosity_plot_min, 18.0)
    log10_viscosity_plot_max = get(kwargs, :log10_viscosity_plot_max, 25.0)
    log10_viscosity_axis_spacing = get(kwargs, :log10_viscosity_axis_spacing, 1.0)
    #figsize = get(kwargs, :figsize, (5.0, 10.0))
    
    fig, ax1 = initialize_plot(figsize_pixels)
    
    # Calculate effective viscosity
    effective_viscosity_no_damage = yield_stress_gridy_no_damage ./ (2.0 * strain_rate)
    effective_viscosity_with_damage = yield_stress_gridy_with_damage ./ (2.0 * strain_rate)
    # Convert to log10
    log10_viscosity_no_damage = log10.(effective_viscosity_no_damage)
    log10_viscosity_with_damage = log10.(effective_viscosity_with_damage)
    # Convert m to km
    depth_km = gridy ./ 1000.0
    
    CairoMakie.lines!(ax1, log10_viscosity_no_damage, depth_km, 
           label="Effective Viscosity (No Damage)", color=:black)
    CairoMakie.lines!(ax1, log10_viscosity_with_damage, depth_km, 
           label="Effective viscosity (Damaged)", color=:black, linestyle=:dot)
    
    ax1.xlabel = "Effective Viscosity log10(Pa.s)"
    ax1.ylabel = "Depth (km)"
    CairoMakie.ylims!(ax1, 0, depth_plot_limit)
    ax1.yticks = get_ticks(depth_plot_limit, depth_axis_spacing)
    ax1.yreversed = true
    CairoMakie.xlims!(ax1, log10_viscosity_plot_min, log10_viscosity_plot_max)

    ax1.xticks = get_ticks(log10_viscosity_plot_max, log10_viscosity_axis_spacing, 
                           log10_viscosity_plot_min)
    
    CairoMakie.axislegend(ax1, position=:rt, framevisible=false, labelsize=10)
    ax1.xgridvisible = true
    ax1.xgridstyle = :dash
    ax1.xgridcolor = :gray
    ax1.xgridwidth = 0.5
    
    text_str = (
        "Strain Rate: $strain_rate 1/s\n" *
        "Cohesion (No Damage): $(cohesions[1]/1e6) MPa\n" *
        "Cohesion (Damaged): $(cohesions[2]/1e6) MPa\n" *
        "Friction Angle (No Damage): $(friction_angles[1]) degrees\n" *
        "Friction Angle (Damaged): $(friction_angles[2]) degrees"
    )
    
    CairoMakie.text!(ax1, text_str, position=(0.45, 0.02), 
          fontsize=7, align=(:left, :bottom))
    
    plot_name = joinpath(plot_output_path, "effective_viscosity.png")
    CairoMakie.save(plot_name, fig)
    
    return nothing
end

""" Initialize plot.

# Arguments
- `figsize_pixels::Tuple{Float64, Float64}`: Figure size
# Returns
- `fig::CairoMakie.Figure`: Figure object
- `ax::CairoMakie.Axis`: Axis object
"""
function initialize_plot(figsize_pixels::Tuple{Float64, Float64})::Tuple{CairoMakie.Figure, CairoMakie.Axis}
    println("Initializing plot with size: $figsize_pixels")
    fig = CairoMakie.Figure(size=figsize_pixels)
    ax = CairoMakie.Axis(fig[1, 1])
    return fig, ax
end

""" Get ticks for axis.

# Arguments
- `scalar_limit::Float64`: Maximum value for ticks
- `spacing::Float64`: Spacing between ticks
- `scalar_min::Float64`: Minimum value for ticks (default: 0.0)
# Returns
- `ticks::Vector{Float64}`: Array of tick values
"""
function get_ticks(
    scalar_limit::Float64, 
    spacing::Float64, 
    scalar_min::Float64=0.0
)::Vector{Float64}
    return collect(scalar_min:spacing:scalar_limit)
end

""" Check if directory exists and create if necessary.

# Arguments
- `path::String`: Directory path to check/create
"""
function check_direc(path::String)::Nothing
    if !isdir(path)
        mkpath(path)
    end
    return nothing
end

end # module
