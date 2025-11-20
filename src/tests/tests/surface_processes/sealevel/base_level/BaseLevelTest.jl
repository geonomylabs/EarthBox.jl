module BaseLevelTest

using Plots
import EarthBox.SurfaceProcesses.Sealevel.RelativeBaseLevel.DensityProps: 
    DensityModel, DensityProperties
import EarthBox.SurfaceProcesses.Sealevel.RelativeBaseLevel: 
    calculate_isostatic_relative_base_level
import EarthBox.SurfaceProcesses.Sealevel.RelativeBaseLevel.ReferenceLithosphere: 
    LithosphereThicknesses, LithosphereThermalProps, reference_lithosphere
import EarthBox.SurfaceProcesses.Sealevel.RelativeBaseLevel.LithostaticPressure: 
    calculate_lithostatic_pressure_from_marker_swarm
import EarthBox.MathTools: linear_interp

struct ComputationalGrid
    xsize::Float64  # x-size of the model in meters
    ysize::Float64  # y-size of the model in meters
    xnum::Int  # x-number of grid nodes
    ynum::Int  # y-number of grid nodes
    sticky_thickness::Float64  # sticky_thickness at top of model domain (meters)
    markers_per_cell_x::Int  # x-spacing of markers in meters
    markers_per_cell_y::Int  # y-spacing of markers in meters
end

function run_test()::Nothing
    # Model grid and swarm parameters
    comp_grid = ComputationalGrid(
        140_000.0,  # xsize
        140_000.0,  # ysize
        160,        # xnum
        160,        # ynum
        10_000.0,   # sticky_thickness
        8,          # markers_per_cell_x
        8           # markers_per_cell_y
    )

    # Analytical solution is calculated only if this is set to true
    turn_off_density_changes = true

    density_model = get_density_model_simple(turn_off_density_changes)

    model_lith_thickness = LithosphereThicknesses(
        (comp_grid.ysize - comp_grid.sticky_thickness),
        22_000.0,
        10_000.0,
        122_000.0,
        100.0
    )

    # If true, use a temperature profile with four linear segments. If false,
    # use the three-layer analytical model that includes a steady-state
    # geotherm in the crust and mantle based on input thermal conductivities
    # and heat production values.
    iuse_linear_segments = 0

    model_lith_thermal = LithosphereThermalProps(
        0.0,
        600.0,
        1330.0,
        0.4,
        2.5,
        2.5,
        2.25,
        9e-7,
        9e-7,
        0.0,
        125_000.0
    )

    (
        marker_x, marker_y, marker_dx, marker_dy, marker_rho,
        gridy, temp_gridy, density_gridy, pressure_gridy
    ) = make_marker_swarm(
        density_model,
        comp_grid,
        model_lith_thickness,
        model_lith_thermal,
        iuse_linear_segments=iuse_linear_segments
    )

    println(
        "Lithostatic pressure at the base of model from grid (GPa): ",
        pressure_gridy[end]/1e9
    )

    y_topo = comp_grid.sticky_thickness
    (
        gridy_for_marker_averaging,
        density_gridy_from_markers,
        pressure_gridy_from_markers
    ) = calculate_lithostatic_pressure_from_marker_swarm(
        marker_x,
        marker_y,
        marker_rho,
        comp_grid.ysize,
        marker_dx*8,
        y_topo,
        marker_dy*4,
        marker_dx*4
    )

    pressure_at_base_of_model_column = pressure_gridy_from_markers[end]

    println(
        "Lithostatic_pressure at the base of model from marker density (GPa): ",
        pressure_at_base_of_model_column/1e9
    )

    reference_lith_thickness = LithosphereThicknesses(
        comp_grid.ysize,
        22_000.0,
        10_000.0,
        122_000.0,
        10.0
    )

    reference_lith_thermal = LithosphereThermalProps(
        0.0,
        600.0,
        1330.0,
        0.4,
        2.5,
        2.5,
        2.25,
        9e-7,
        9e-7,
        0.0,
        125_000.0
    )

    (
        relative_base_level,
        gridy_ref,
        temp_gridy_ref,
        density_gridy_ref,
        pressure_gridy_ref
    ) = calculate_isostatic_relative_base_level(
        model_lith_thickness.total_column_thickness_meters,
        pressure_at_base_of_model_column,
        reference_lith_thickness,
        reference_lith_thermal,
        density_model,
        iuse_linear_segments=iuse_linear_segments
    )

    println(
        "Lithostatic pressure at the base of reference lithosphere (GPa): ",
        pressure_gridy_ref[end]/1e9
    )

    println("relative_base_level (m): ", relative_base_level)

    thickness_cont_crust_meters = (
        model_lith_thickness.thickness_upr_cont_crust_meters
        + model_lith_thickness.thickness_lwr_cont_crust_meters
    )
    thickness_cont_crust_ref_meters = (
        reference_lith_thickness.thickness_upr_cont_crust_meters
        + reference_lith_thickness.thickness_lwr_cont_crust_meters
    )

    if turn_off_density_changes
        relative_base_level_analytical = analytical_solution(
            density_model,
            thickness_cont_crust_meters,
            thickness_cont_crust_ref_meters
        )

        println(
            "relative_base_level_analytical (m): ",
            relative_base_level_analytical
        )
    end

    make_temperature_plots(
        gridy, temp_gridy,
        gridy_ref, temp_gridy_ref
    )
    make_density_plots(
        gridy, density_gridy,
        gridy_ref, density_gridy_ref,
        gridy_for_marker_averaging, density_gridy_from_markers,
        show_plot=false
    )
    make_pressure_plots(
        gridy, pressure_gridy,
        gridy_ref, pressure_gridy_ref,
        gridy_for_marker_averaging, pressure_gridy_from_markers,
        show_plot=false
    )

    return nothing
end

function make_marker_swarm(
    density_model::DensityModel,
    comp_grid::ComputationalGrid,
    model_lith_thickness::LithosphereThicknesses,
    model_lith_thermal::LithosphereThermalProps;
    iuse_linear_segments::Int64=1
)::Tuple{Vector{Float64}, Vector{Float64}, Float64, Float64, Vector{Float64},
         Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}
    
    (
        marker_x, marker_y, marker_dx, marker_dy
    ) = calculate_marker_coordinates(
        comp_grid.xnum,
        comp_grid.ynum,
        comp_grid.markers_per_cell_x,
        comp_grid.markers_per_cell_y,
        comp_grid.xsize,
        comp_grid.ysize
    )

    (
        gridy, temp_gridy, density_gridy, pressure_gridy
    ) = reference_lithosphere(
        density_model,
        model_lith_thickness,
        model_lith_thermal,
        number_of_pressure_iterations=3,
        iuse_linear_segments=iuse_linear_segments
    )

    marker_rho = calculate_marker_density(
        gridy, density_gridy, marker_y, comp_grid.sticky_thickness
    )

    return (
        marker_x, marker_y, marker_dx, marker_dy, marker_rho,
        gridy, temp_gridy, density_gridy, pressure_gridy
    )
end

function calculate_marker_coordinates(
    xnum::Int,
    ynum::Int,
    markers_per_cell_x::Int,
    markers_per_cell_y::Int,
    xsize::Float64,
    ysize::Float64
)::Tuple{Vector{Float64}, Vector{Float64}, Float64, Float64}
    
    n_markers_x = xnum * markers_per_cell_x
    n_markers_y = ynum * markers_per_cell_y

    nmarkers = n_markers_x * n_markers_y

    marker_x = zeros(nmarkers)
    marker_y = zeros(nmarkers)

    marker_dx = xsize / xnum / markers_per_cell_x
    marker_dy = ysize / ynum / markers_per_cell_y

    marker_counter = 1
    for j in 1:n_markers_y
        for i in 1:n_markers_x
            marker_x[marker_counter] = (i - 1) * marker_dx
            marker_y[marker_counter] = (j - 1) * marker_dy
            marker_counter += 1
        end
    end

    return marker_x, marker_y, marker_dx, marker_dy
end

function calculate_marker_density(
    gridy::Vector{Float64},
    density_gridy::Vector{Float64},
    marker_y::Vector{Float64},
    sticky_thickness::Float64
)::Vector{Float64}
    
    nnodes = length(gridy)
    min_y = gridy[1]
    max_y = gridy[end]
    y_index_max = nnodes

    nmarkers = length(marker_y)
    marker_rho = zeros(nmarkers)
    for imarker in 1:nmarkers
        y = marker_y[imarker] - sticky_thickness
        marker_rho[imarker] = linear_interp(
            y, gridy, density_gridy, nnodes, min_y, max_y, y_index_max
        )
    end
    return marker_rho
end

function make_temperature_plots(
    gridy::Vector{Float64},
    temp_gridy::Vector{Float64},
    gridy_ref::Vector{Float64},
    temp_gridy_ref::Vector{Float64};
    show_plot::Bool=false
)::Nothing
    
    p = plot(
        temp_gridy_ref .- 273, gridy_ref./1000.0,
        label="Temperature_Ref",
        color=:red
    )
    plot!(
        temp_gridy .- 273, gridy./1000.0,
        label="Temperature_Model_Grid",
        color=:blue
    )
    xlabel!("Temperature (C)")
    ylabel!("y (km below rock line)")
    xlims!(0, 1400)
    ylims!(0, 160)
    yaxis!(:flip)
    
    if !show_plot
        savefig(p, "temperature_profiles.png")
    else
        display(p)
    end
    
    return nothing
end

function make_density_plots(
    gridy::Vector{Float64},
    density_gridy::Vector{Float64},
    gridy_ref::Vector{Float64},
    density_gridy_ref::Vector{Float64},
    gridy_for_marker_averaging::Vector{Float64},
    density_gridy_from_markers::Vector{Float64};
    show_plot::Bool=false
)::Nothing
    
    p = plot(
        density_gridy_ref, gridy_ref./1000.0,
        label="Density_Ref",
        color=:red,
        legend=:bottomleft
    )
    plot!(
        density_gridy, gridy./1000.0,
        label="Density_Model_Grid",
        color=:blue
    )
    plot!(
        density_gridy_from_markers, gridy_for_marker_averaging./1000.0,
        label="Density_Marker_Interpolation_Grid",
        color=:green
    )
    xlabel!("Density (kg/m3)")
    ylabel!("y (km below rock line)")
    xlims!(2500, 3500)
    ylims!(0, 160)
    yaxis!(:flip)
    
    if !show_plot
        savefig(p, "density_profiles.png")
    else
        display(p)
    end
    
    return nothing
end

function make_pressure_plots(
    gridy::Vector{Float64},
    pressure_gridy::Vector{Float64},
    gridy_ref::Vector{Float64},
    pressure_gridy_ref::Vector{Float64},
    gridy_for_marker_averaging::Vector{Float64},
    pressure_gridy_from_markers::Vector{Float64};
    show_plot::Bool=false
)::Nothing
    
    p = plot(
        pressure_gridy_ref./1e9, gridy_ref./1000.0,
        label="Pressure_Ref",
        color=:red
    )
    plot!(
        pressure_gridy./1e9, gridy./1000.0,
        label="Pressure_Model_Grid",
        color=:blue
    )
    plot!(
        pressure_gridy_from_markers./1e9, gridy_for_marker_averaging./1000.0,
        label="Pressure_Marker_Interpolation_Grid",
        color=:green
    )
    xlabel!("Pressure (GPa)")
    ylabel!("y (km below rock line)")
    xlims!(0, 6.0)
    ylims!(0, 160)
    yaxis!(:flip)
    
    if !show_plot
        savefig(p, "pressure_profiles.png")
    else
        display(p)
    end
    
    return nothing
end

function analytical_solution(
    density_model::DensityModel,
    thickness_cont_crust_meters::Float64,
    thickness_cont_crust_reference_meters::Float64;
    density_water::Float64=1000.0
)::Float64
    
    density_crust = get_density(density_model.upr_cont_crust)
    density_mantle = get_density(density_model.asthenosphere)

    relative_base_level_analytical = (
        (thickness_cont_crust_reference_meters - thickness_cont_crust_meters)
        * (density_mantle - density_crust)
        / (density_mantle - density_water)
    )
    return -relative_base_level_analytical
end

function get_density(density_props::DensityProperties)::Float64
    return density_props.standard_density
end

function get_density_model_simple(
    turn_off_density_change::Bool=true
)::DensityModel
    
    (
        density_props_upr_cont_crust,
        density_props_lwr_cont_crust,
    ) = get_density_properties_continental_crust_simple(turn_off_density_change)

    (
        density_props_upper_mantle_lithosphere,
        density_props_middle_mantle_lithosphere,
        density_props_lower_mantle_lithosphere
    ) = get_density_properties_mantle_lithosphere_simple(turn_off_density_change)

    (
        density_props_asthenosphere
    ) = get_density_properties_asthenosphere_simple(turn_off_density_change)

    density_model = DensityModel(
        density_props_upr_cont_crust,
        density_props_lwr_cont_crust,
        density_props_upper_mantle_lithosphere,
        density_props_middle_mantle_lithosphere,
        density_props_lower_mantle_lithosphere,
        density_props_asthenosphere
    )

    return density_model
end

function get_density_properties_asthenosphere_simple(
    turn_off_density_change::Bool=true
)::DensityProperties
    
    if turn_off_density_change
        expansivity = 0.0
        compressibility = 0.0
    else
        expansivity = 3.0e-5
        compressibility = 1.0e-11
    end

    density_props_asthenosphere = DensityProperties(3300.0, expansivity, compressibility)
    return density_props_asthenosphere
end

function get_density_properties_mantle_lithosphere_simple(
    turn_off_density_change::Bool=true
)::Tuple{DensityProperties, DensityProperties, DensityProperties}
    
    if turn_off_density_change
        expansivity = 0.0
        compressibility = 0.0
    else
        expansivity = 3.0e-5
        compressibility = 1.0e-11
    end

    density_props_upr_mantle_lith = DensityProperties(3300.0, expansivity, compressibility)
    density_props_mid_mantle_lith = DensityProperties(3300.0, expansivity, compressibility)
    density_props_lwr_mantle_lith = DensityProperties(3300.0, expansivity, compressibility)
    return density_props_upr_mantle_lith, density_props_mid_mantle_lith, density_props_lwr_mantle_lith
end

function get_density_properties_continental_crust_simple(
    turn_off_density_change::Bool=true
)::Tuple{DensityProperties, DensityProperties}
    
    if turn_off_density_change
        expansivity = 0.0
        compressibility = 0.0
    else
        expansivity = 3.0e-5
        compressibility = 1.0e-11
    end

    density_props_upr_cont_crust = DensityProperties(2800.0, expansivity, compressibility)
    density_props_lwr_cont_crust = DensityProperties(2800.0, expansivity, compressibility)
    return density_props_upr_cont_crust, density_props_lwr_cont_crust
end

end # module 