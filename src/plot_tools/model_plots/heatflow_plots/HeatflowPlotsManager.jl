module HeatflowPlotsManager

using Printf
import CairoMakie
import EarthBox.Markers.MarkerMaterials.MaterialsContainer: Materials
import EarthBox.Markers.MarkerMaterials.GetMaterialIDs: get_sediment_material_id
import EarthBox.Markers.MarkerMaterials.GetMaterialIDs: get_solidified_basalt_material_id
import EarthBox.JLDTools: intstr
import EarthBox.JLDTools: get_jld_data, get_jld_topo_data
import EarthBox.JLDTools: get_jld_marker_id_data
import EarthBox.JLDTools: get_marker_coordinate_arrays
import EarthBox.MathTools: linear_interp_vals_opt!
import EarthBox.MathTools: linear_interp_bisection
import EarthBox.MathTools: smooth_surface
import EarthBox.ModelStructureManager.TopAndBottom: calculate_top_and_bottom_of_layer_opt
import ..PlotDtypes: AxesType
import ..PlotTimeSteppingManager: PlotTimeStepping
import ..GridPlotsManager.DataNames: ScalarH5Datanames
#import ..MarkerPlotsManager.DataNames: MarkerDataNames # Circular import, need to refactor code

"""
    HeatflowPlots

Mutable struct to manage heat flow plots and associated data.

# Fields
- `output_dir_path::String`: Path to output directory
- `plot_output_path::String`: Path for plot outputs
- `scalar_datanames::ScalarH5Datanames`: Scalar H5 (now JLD) data names
- `time_stepping::PlotTimeStepping`: Plot time stepping configuration
- `ioutput::Union{Int, Nothing}`: Current output index
- `marker_data_names::Union{MarkerDataNames, Nothing}`: Marker data names
"""
mutable struct HeatflowPlots
    output_dir_path::String
    plot_output_path::String
    scalar_datanames::ScalarH5Datanames
    time_stepping::PlotTimeStepping
    ioutput::Union{Int, Nothing}
    marker_data_names::Union{Any, Nothing} # Change Any to MarkerDataNames after refactoring code
end

function HeatflowPlots(
    output_dir_path::String,
    plot_output_path::String,
    time_stepping::PlotTimeStepping,
    marker_data_names::Union{Any, Nothing} = nothing # Change Any to MarkerDataNames after refactoring code
)
    scalar_datanames = ScalarH5Datanames()
    return HeatflowPlots(
        output_dir_path, plot_output_path, scalar_datanames, 
        time_stepping, nothing, marker_data_names
    )
end

function set_ioutput!(manager::HeatflowPlots, ioutput::Int)::Nothing
    manager.ioutput = ioutput
    return nothing
end

function get_heat_flow_grids(
    manager::HeatflowPlots;
    materials::Union{Materials, Nothing} = nothing,
    use_sediment_thickness::Bool = false
)::Tuple{
    Vector{Float64}, Vector{Float64}, Vector{Float64}, 
    Vector{Float64}, Vector{Float64}, Vector{Float64}
    }
    # Get temperature grid
    (_model_time, temp_gridy, temp_gridx, temp_array) = get_temperature_grid(manager)
    
    # Get thermal conductivity grid
    (_model_time, _therm_cond_gridy, _therm_cond_gridx, therm_cond_array) = 
        get_therm_cond_grid(manager)
    
    # Get topography data
    (topo_gridy, topo_gridx, _length_units) = get_topo_data(manager)
    
    # Calculate overburden thickness
    if !use_sediment_thickness || materials === nothing
        overburden_thickness_topo_gridx = zeros(length(topo_gridx))
    else
        overburden_thickness_topo_gridx = calculate_sediment_thickness(
            manager, topo_gridx, materials)
    end
    
    # Calculate heat flow
    (heat_flow_x, temperature1_x, temperature2_x, conductivity1_x, 
     conductivity2_x) = calculate_heat_flow(
        temp_gridx, temp_gridy, therm_cond_array, temp_array,
        topo_gridx, topo_gridy, overburden_thickness_topo_gridx
    )
    
    # Smooth the heat flow surface
    heat_flow_x = smooth_surface(heat_flow_x, nsmooth=8)
    
    return (heat_flow_x, temp_gridx, temperature1_x, temperature2_x, 
            conductivity1_x, conductivity2_x)
end

function get_temperature_grid(
    manager::HeatflowPlots
)::Tuple{Float64, Vector{Float64}, Vector{Float64}, Matrix{Float64}}
    jld_filename = get_jld_field_filename(manager)
    dataname = manager.scalar_datanames.TempC
    (
        model_time, temp_gridy, temp_gridx, temp_array, _units_dict
    ) = get_jld_data(dataname, jld_filename)
    return model_time, temp_gridy, temp_gridx, temp_array
end

function get_therm_cond_grid(
    manager::HeatflowPlots;
    conductivity_max::Float64 = 2.5
)::Tuple{Float64, Vector{Float64}, Vector{Float64}, Matrix{Float64}}
    jld_filename = get_jld_field_filename(manager)
    dataname = manager.scalar_datanames.therm_cond
    (model_time, therm_cond_gridy, therm_cond_gridx, therm_cond_array, 
     _units_dict) = get_jld_data(dataname, jld_filename)
    # Cap conductivity values
    therm_cond_array = min.(therm_cond_array, conductivity_max)
    return model_time, therm_cond_gridy, therm_cond_gridx, therm_cond_array
end

function get_jld_field_filename(manager::HeatflowPlots)::String
    output_dir_path = manager.output_dir_path
    ioutput = manager.ioutput
    jld_filename = joinpath(output_dir_path, "fields_$(intstr(ioutput)).jld")
    return jld_filename
end

function get_topo_data(
    manager::HeatflowPlots
)::Tuple{Vector{Float64}, Vector{Float64}, String}
    jld_filename = get_jld_topo_filename(manager)
    (topo_gridy, topo_gridx, length_units) = get_jld_topo_data(jld_filename)
    return topo_gridy, topo_gridx, length_units
end

function get_jld_topo_filename(manager::HeatflowPlots)::String
    output_dir_path = manager.output_dir_path
    ioutput = manager.ioutput
    jld_filename = joinpath(output_dir_path, "topo_$(intstr(ioutput)).jld")
    return jld_filename
end

function calculate_sediment_thickness(
    manager::HeatflowPlots,
    topo_gridx::Vector{Float64},
    materials::Materials
)::Vector{Float64}
    (
        _model_time_myr, marker_x_m, marker_y_m, marker_matids
    ) = get_jld_marker_data(manager)
    
    search_radius = 0.5 * (topo_gridx[2] - topo_gridx[1])
    
    matid_sed = get_sediment_material_id(materials)
    matid_basalt = get_solidified_basalt_material_id(materials)
    
    material_ids_of_layer = Int16[-1, -1]
    material_ids_of_layer[1] = matid_sed
    material_ids_of_layer[2] = matid_basalt
    
    (top, bottom) = calculate_top_and_bottom_of_layer_opt(
        material_ids_of_layer, marker_matids, marker_x_m, marker_y_m,
        topo_gridx, search_radius; use_smoothing=true, nsmooth=2
    )
    
    sediment_thickness_topo_gridx = (bottom .- top) ./ 1000.0
    
    return sediment_thickness_topo_gridx
end

function get_jld_marker_data(
    manager::HeatflowPlots
)::Tuple{Float64, Vector{Float64}, Vector{Float64}, Vector{Int16}}
    file_base_name = "particles"
    (
        marker_x_m, marker_y_m, marker_matid, model_time_myr
    ) = read_marker_material_ids(manager.output_dir_path, file_base_name, manager.ioutput)
    return model_time_myr, marker_x_m, marker_y_m, marker_matid
end

# This function is duplicated from the benchmarks utils module Reader.jl
function read_marker_material_ids(
    input_path::String,
    _base_name::String,
    ioutput::Int
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Int16}, Float64}
    tmyr, marker_matid = get_jld_marker_id_data(ioutput, input_path)
    tmyr, marker_x_m, marker_y_m = get_marker_coordinate_arrays(ioutput, input_path)
    return marker_x_m, marker_y_m, marker_matid, tmyr
end

function plot_heatflow(
    axes::CairoMakie.Axis,
    gridx::Vector{Float64},
    heat_flow::Vector{Float64},
    heat_flow_basal::Vector{Float64};
    iplot_basal_heat_flow::Bool = false,
    legendfontsize::Int = 10
)::Nothing
    CairoMakie.lines!(axes, gridx, heat_flow; label="Surface Heat Flow", color="red")
    if iplot_basal_heat_flow
        CairoMakie.lines!(axes, gridx, heat_flow_basal; 
              label="Basal Heat Flow", color="blue", linestyle=:dash)
    end
    axes.ylabel = "Heat Flow (mW/m²)"
    CairoMakie.axislegend(
        axes, position=:rt, labelsize=legendfontsize, patchsize=(5,5), rowgap=2)
    return nothing
end

function plot_heatflow_tests(manager::HeatflowPlots)
    for ioutput in manager.time_stepping.steps
        make_test_plots(manager, ioutput)
    end
end

function make_test_plots(
    manager::HeatflowPlots,
    ioutput::Int;
    min_hf::Float64 = 20.0,
    max_hf::Float64 = 200.0,
    hf_spacing::Float64 = 10.0,
    plot_conductivities::Bool = false,
    plot_temperatures::Bool = false
)
    set_ioutput!(manager, ioutput)
    
    (heat_flow_x, temp_gridx, temperature1_x, temperature2_x,
     conductivity1_x, conductivity2_x) = get_heat_flow_grids(manager)
    
    # Main heat flow plot
    fig = CairoMakie.Figure(size=(1000, 500))
    ax = CairoMakie.Axis(fig[1, 1])
    CairoMakie.lines!(ax, temp_gridx, heat_flow_x; color="black")
    CairoMakie.xlims!(ax, 0, temp_gridx[end])
    CairoMakie.ylims!(ax, min_hf, max_hf)
    CairoMakie.xlabel!(ax, "Distance (km)")
    CairoMakie.ylabel!(ax, "Heat Flow (mW/m)")
    CairoMakie.title!(ax, "Heat Flow")
    
    CairoMakie.save(joinpath(manager.plot_output_path, "heat_flow_$(manager.ioutput).png"), 
         fig)
    
    if plot_temperatures
        fig_temp = CairoMakie.Figure()
        ax_temp = CairoMakie.Axis(fig_temp[1, 1])
        CairoMakie.lines!(ax_temp, temp_gridx, temperature1_x; color="black")
        CairoMakie.lines!(ax_temp, temp_gridx, temperature2_x; color="red")
        CairoMakie.xlims!(ax_temp, 0, temp_gridx[end])
        CairoMakie.ylims!(ax_temp, 0, 500.0)
        CairoMakie.xlabel!(ax_temp, "Distance (km)")
        CairoMakie.ylabel!(ax_temp, "Temperature (C)")
        CairoMakie.title!(ax_temp, "Temperature")
        
        CairoMakie.save(joinpath(manager.plot_output_path, 
                      "temperature_$(manager.ioutput).png"), fig_temp)
    end
    
    if plot_conductivities
        fig_cond = CairoMakie.Figure()
        ax_cond = CairoMakie.Axis(fig_cond[1, 1])
        CairoMakie.lines!(ax_cond, temp_gridx, conductivity1_x; color="black")
        CairoMakie.lines!(ax_cond, temp_gridx, conductivity2_x; color="red")
        CairoMakie.xlims!(ax_cond, 0, temp_gridx[end])
        CairoMakie.ylims!(ax_cond, 1.0, 3.0)
        CairoMakie.xlabel!(ax_cond, "Distance (km)")
        CairoMakie.ylabel!(ax_cond, "Conductivity (W/mK)")
        CairoMakie.title!(ax_cond, "Conductivity")
        
        CairoMakie.save(joinpath(manager.plot_output_path, 
                      "conductivity_$(manager.ioutput).png"), fig_cond)
    end
end

function calculate_heat_flow(
    gridx::Vector{Float64},
    gridy::Vector{Float64},
    conductivity_grid::Matrix{Float64},
    temperature_grid::Matrix{Float64},
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64},
    overburden_thickness_topo_gridx::Vector{Float64};
    dzz_hf::Float64 = 0.25,
    z_shift_hf::Float64 = 1.0
)
    xnum = length(gridx)
    ynum = length(gridy)
    
    heat_flow_x = zeros(xnum)
    temperature1_x = zeros(xnum)
    temperature2_x = zeros(xnum)
    conductivity1_x = zeros(xnum)
    conductivity2_x = zeros(xnum)
    
    # Note: Julia uses 1-based indexing and column-major order
    for j in 1:xnum
        x_location = gridx[j]
        
        y_topo = linear_interp_bisection(
            topo_gridx ./ 1000.0, topo_gridy ./ 1000.0, x_location)
        
        overburden = linear_interp_bisection(
            topo_gridx ./ 1000.0, overburden_thickness_topo_gridx, x_location)
        
        # Extract temperature and conductivity profiles for this x location
        temperature_y = zeros(ynum)
        conductivity_y = zeros(ynum)
        for i in 1:ynum
            temperature_y[i] = temperature_grid[i, j]
            conductivity_y[i] = conductivity_grid[i, j]
        end
        
        z_shift_hf = (gridy[2] - gridy[1]) * 4.0 + overburden
        dzz_hf = (gridy[2] - gridy[1]) * 2.0
        
        (
            heat_flow_x[j], temperature1_x[j], temperature2_x[j],
            conductivity1_x[j], conductivity2_x[j]
        ) = calc_heat_flow(y_topo, gridy, temperature_y, conductivity_y, dzz_hf, z_shift_hf)
    end
    
    return (heat_flow_x, temperature1_x, temperature2_x, 
            conductivity1_x, conductivity2_x)
end

"""
    calc_heat_flow(y_topo::Float64, gridy::Vector{Float64},
                  temperature_y::Vector{Float64}, conductivity_y::Vector{Float64},
                  dzz::Float64, z_shift::Float64)

Calculate heat flow at a given depth y.

# Arguments
- `y_topo::Float64`: Topography y coordinate
- `gridy::Vector{Float64}`: Grid y coordinates
- `temperature_y::Vector{Float64}`: Temperature profile
- `conductivity_y::Vector{Float64}`: Conductivity profile
- `dzz::Float64`: Depth spacing
- `z_shift::Float64`: Depth shift

# Returns
- `Tuple`: (hf, temperature_1, temperature_2, conductivity_1, conductivity_2)
"""
function calc_heat_flow(
    y_topo::Float64,
    gridy::Vector{Float64},
    temperature_y::Vector{Float64},
    conductivity_y::Vector{Float64},
    dzz::Float64,
    z_shift::Float64
)
    ny_interp = 2
    
    y_interp = zeros(ny_interp)
    y_interp[1] = y_topo + z_shift
    y_interp[2] = y_topo + dzz + z_shift
    
    vals_interp = zeros(ny_interp)
    linear_interp_vals_opt!(gridy, temperature_y, y_interp, vals_interp)
    
    temperature_1 = vals_interp[1]
    temperature_2 = vals_interp[2]
    
    vals_interp = zeros(ny_interp)
    linear_interp_vals_opt!(gridy, conductivity_y, y_interp, vals_interp)
    
    conductivity_1 = vals_interp[1]
    conductivity_2 = vals_interp[2]
    
    conductivity = (conductivity_1 + conductivity_2) / 2.0
    
    # Convert to mW/m²
    hf = conductivity * (temperature_2 - temperature_1) / (dzz * 1000.0) * 1000.0
    
    return hf, temperature_1, temperature_2, conductivity_1, conductivity_2
end

end # module
