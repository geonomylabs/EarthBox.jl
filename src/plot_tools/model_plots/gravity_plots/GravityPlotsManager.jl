module GravityPlotsManager

using Printf
import CairoMakie
import JLD2
import EarthBox.JLDTools: intstr
import EarthBox.JLDTools: get_jld_data, get_jld_topo_data
import EarthBox.JLDTools: get_y_sealevel_from_jld_marker_data
import EarthBox.JLDTools: get_base_level_shift_from_jld_marker_data
import EarthBox.Gravity: calculate_free_air_and_bouguer
import EarthBox.MathTools: smooth_surface
import ..PlotDtypes: AxesType
import ..PlotTimeSteppingManager: PlotTimeStepping
import ..GridPlotsManager.DataNames: ScalarH5Datanames

mutable struct GravityPlots
    output_dir_path::String
    plot_output_path::String
    scalar_datanames::ScalarH5Datanames
    time_stepping::PlotTimeStepping
    ioutput::Union{Int, Nothing}
end

function GravityPlots(
    output_dir_path::String,
    plot_output_path::String,
    time_stepping::PlotTimeStepping
)
    scalar_datanames = ScalarH5Datanames()
    return GravityPlots(
        output_dir_path, plot_output_path, scalar_datanames, 
        time_stepping, nothing
    )
end

function set_ioutput!(manager::GravityPlots, ioutput::Int)::Nothing
    manager.ioutput = ioutput
    return nothing
end

function get_gravity_grids(
    manager::GravityPlots
)::Tuple{Float64, Vector{Float64}, Vector{Float64}, Vector{Float64}}
    """Calculate gravity grid at model gridx nodes and associated grids.

    Return
    ------
    model_time : Float64
        Model time (Myr).

    topo_gridx : Vector{Float64}
        Topography x (m).

    gravity_grid_mgal : Vector{Float64}
        Gravity grid (mgal).

    gravity_grid_free_air_mgal : Vector{Float64}
        Gravity free air grid (mgal).
    """
    (model_time, gridy, gridx, rho_array) = get_density_grid(manager)
    (topo_gridy, topo_gridx, _length_units) = get_topo_data(manager)
    iout_tag = manager.ioutput === nothing ? "none" : intstr(manager.ioutput)

    println(">> Density min: ", minimum(rho_array))
    println(">> Density max: ", maximum(rho_array))

    println(">> gridx min: ", minimum(gridx))
    println(">> gridx max: ", maximum(gridx))

    println(">> gridy min: ", minimum(gridy))
    println(">> gridy max: ", maximum(gridy))

    println(">> topo_gridx min: ", minimum(topo_gridx))
    println(">> topo_gridx max: ", maximum(topo_gridx))

    y_sealevel = get_y_sealevel_from_jld_marker_data(
        manager.ioutput, manager.output_dir_path)

    println(">> y_sealevel: ", y_sealevel)

    base_level_shift = get_base_level_shift_from_jld_marker_data(
        manager.ioutput, manager.output_dir_path)
    println(">> Base level shift from jld data: ", base_level_shift)
    # Remove the base level shift to get to the global sea level
    y_sealevel = y_sealevel - base_level_shift

    println(">> y_sealevel (global): ", y_sealevel)

    (gravity_grid_mgal, gravity_grid_free_air_mgal) = 
        calculate_free_air_and_bouguer(
            gridx, gridy, rho_array, y_sealevel, topo_gridx, topo_gridy)

    println(">> Gravity free air/bouguer min: ", minimum(gravity_grid_mgal))
    println(">> Gravity free air/bouguer max: ", maximum(gravity_grid_mgal))

    println(">> Gravity free air min: ", minimum(gravity_grid_free_air_mgal))
    println(">> Gravity free air max: ", maximum(gravity_grid_free_air_mgal))

    return model_time, topo_gridx, gravity_grid_mgal, gravity_grid_free_air_mgal
end


function make_jld2_gravity_file(
    gridx_km::Vector{Float64}, 
    gravity_grid_mgal::Vector{Float64}, 
    gravity_grid_free_air_mgal::Vector{Float64},
    output_dir::String,
    ioutput::Int,
    model_time::Float64,
    time_units::String
)::Nothing
    jld_gravity_filename = get_jld_gravity_filename(output_dir, ioutput)
    jld_gravity_file_path = joinpath(output_dir, jld_gravity_filename)
    
    JLD2.jldopen(jld_gravity_file_path, "w") do file
        file["Myr"] = model_time
        file["noutput"] = ioutput
        
        group = JLD2.Group(file, "gridx")
        group["array"] = gridx_km
        group["name"] = "gridx"
        group["units"] = "km"
        group["noutput"] = ioutput
        group["model_time"] = model_time
        group["time_units"] = time_units
        
        group = JLD2.Group(file, "gravity_bouguer_mgal")
        group["array"] = gravity_grid_mgal
        group["name"] = "gravity_bouguer_mgal"
        group["units"] = "mgal"
        group["noutput"] = ioutput
        group["model_time"] = model_time
        group["time_units"] = time_units

        group = JLD2.Group(file, "gravity_free_air_mgal")
        group["array"] = gravity_grid_free_air_mgal
        group["name"] = "gravity_free_air_mgal"
        group["units"] = "mgal"
        group["noutput"] = ioutput
        group["model_time"] = model_time
        group["time_units"] = time_units

    end
    return nothing
end

function get_jld_gravity_filename(output_dir_path::String, ioutput::Int)::String
    """Get the jld gravity filename."""
    jld_filename = joinpath(output_dir_path, "gravity_$(intstr(ioutput)).jld")
    return jld_filename
end

function get_density_grid(
    manager::GravityPlots
)::Tuple{Float64, Vector{Float64}, Vector{Float64}, Matrix{Float64}}
    """Get density grid.

    Returns
    -------
    model_time : Float64
        Model time (Myr).

    gridy : Vector{Float64}
        Density grid y (km).

    gridx : Vector{Float64}
        Density grid x (km).

    rho_array : Matrix{Float64}
        Density array (kg/m³).
    """
    jld_filename = get_jld_field_filename(manager)

    dataname = manager.scalar_datanames.rho_kg_m3
    (model_time, gridy, gridx, rho_array, _units_dict) = 
        get_jld_data(dataname, jld_filename)
    gridx = gridx .* 1000.0
    gridy = gridy .* 1000.0
    return model_time, gridy, gridx, rho_array
end

function get_jld_field_filename(manager::GravityPlots)::String
    """Get jld field filename."""
    output_dir_path = manager.output_dir_path
    ioutput = manager.ioutput
    jld_filename = joinpath(output_dir_path, "fields_$(intstr(ioutput)).jld")
    return jld_filename
end

function get_topo_data(
    manager::GravityPlots
)::Tuple{Vector{Float64}, Vector{Float64}, String}
    """Get topography data.

    Returns
    -------
    topoy : Vector{Float64}
        Topography y (m).

    topox : Vector{Float64}
        Topography x (m).

    length_units : String
        Length units.
    """
    jld_filename = get_jld_topo_filename(manager)
    (topo_gridy, topo_gridx, length_units) = get_jld_topo_data(jld_filename)
    return topo_gridy, topo_gridx, length_units
end

function get_jld_topo_filename(manager::GravityPlots)::String
    """Get the jld topography filename."""
    output_dir_path = manager.output_dir_path
    ioutput = manager.ioutput
    jld_filename = joinpath(output_dir_path, "topo_$(intstr(ioutput)).jld")
    return jld_filename
end


function plot_gravity(
    axes::CairoMakie.Axis,
    gridx::Vector{Float64},
    gravity_grid_mgal::Vector{Float64},
    gravity_grid_free_air_mgal::Vector{Float64};
    legendfontsize::Int = 10
)::Nothing
    """Plot gravity."""
    CairoMakie.lines!(axes, gridx, gravity_grid_mgal;
        label="Free Air (Offshore)/ Bouguer (Onshore)", color="red")
    CairoMakie.lines!(axes, gridx, gravity_grid_free_air_mgal;
        label="Free Air", color="blue", linestyle=:dash)
    axes.ylabel = "Gravity Anomaly (mgal)"
    CairoMakie.axislegend(
        axes, position=:rt, labelsize=legendfontsize, patchsize=(5,5), rowgap=2)
    return nothing
end

function plot_density_test(manager::GravityPlots, rho_grid::Matrix{Float64})
    """Plot density grid for testing."""
    fig = CairoMakie.Figure(size=(1000, 800))
    ax = CairoMakie.Axis(fig[1, 1])
    
    # Note: Julia uses column-major order, so we need to transpose for display
    hm_test = CairoMakie.heatmap!(ax, rho_grid'; colormap=:viridis)
    CairoMakie.Colorbar(fig[1, 2], hm_test; label="Density")
    ax.title = "2D Density Grid"
    ax.xlabel = "X"
    ax.ylabel = "Y"
    
    CairoMakie.save(joinpath(manager.plot_output_path, "density_grid.png"), 
        fig; px_per_unit=3)
    return nothing
end

end # module
