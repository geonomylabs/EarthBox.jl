module SubgridDiffusion

import EarthBox.ModelDataContainer: ModelData
import EarthBox: GridFuncs
import EarthBox.Interpolation: WeightFuncs
import EarthBox.Interpolation: BilinearAverage
import EarthBox.Interpolation: GridToMarker

""" Compute and apply subgrid changes for markers and interpolate to grid.

This function calculates subgrid temperature changes for markers delta_temp_subgrid, 
updates the marker temperature array `marker_TK` by adding the subgrid temperature 
change, and then interpolates the subgrid temperature change to the basic grid 
array `dtkn`.

Marker temperature is updated using the following equation:

    T_m = T_m + dT_sg,                                                 eq. 1

where T_m is the current marker temperature and dT_sg is the subgrid temperature 
change defined as:

    dT_sg = dT_sg_tot*(1-exp(C_sg)),                                   eq. 2

and dT_sg_tot is the total subgrid temperature change defined as:

    dT_sg_tot = (T_nodal_m - T_m),                                     eq. 3

where T_nodal_m is the current grid temperature array tk1 interpolated to the 
marker. The term f_sg is defined as follows:

    f_sg = -D*dt/t_sg,                                                 eq. 4

where D is the subgrid diffusion coefficient, dt is the current timestep, and 
t_sg is the thermal diffusion timescale defined as:
    
    t_sg = rho_m*Cp_m/k_m/(2.0/xstpavg^2.0 + 2.0/ystpavg^2.0)        eq. 5

where rho_m is the density of the marker, Cp_m is the heat capacity of the marker, 
k is the conductivity of the marker, and xstpavg and ystpavg are the average 
x- and y-spacing of the basic grid.

The parameters T_nodal_m, k_m and rho_m*Cp_m are interpolated to the marker from 
the basic grid arrays `tk1`, `kt1` and `rhocp1`, respectively, using bilinear 
interpolation.

# Arguments
- `model::ModelData`: The model data container

# Updated Arrays
## model.markers.arrays.thermal
- marker_TK.array: Array{Float64,1}
    - Marker temperature (K).

## model.heat_equation.arrays.temperature
- dtkn.array: Array{Float64,2}
    - Subgrid temperature change on basic grid interpolated from markers.
"""
function subgrid_heat_diffusion_update!(model::ModelData, inside_flags::Vector{Int8})::Nothing
    (
        delta_temp_subgrid_markers
    ) = calculate_subgrid_temperature_change_and_correct_marker_temperature!(model, inside_flags)
    interpolate_subgrid_temperature_change_to_grid!(model, delta_temp_subgrid_markers, inside_flags)
    return nothing
end

function calculate_subgrid_temperature_change_and_correct_marker_temperature!(
    model::ModelData, 
    inside_flags::Vector{Int8},
)::Vector{Float64}
    marknum::Int64 = model.markers.parameters.distribution.marknum.value
    xstpavg::Float64 = model.grids.parameters.geometry.xstpavg.value
    ystpavg::Float64 = model.grids.parameters.geometry.ystpavg.value
    timestep::Float64 = model.timestep.parameters.main_time_loop.timestep.value
    diffusion_coef::Float64 = model.markers.parameters.subgrid_diffusion.subgrid_diff_coef_temp.value

    marker_xn::Vector{Int32} = model.markers.arrays.grid_marker_relationship.marker_xn.array
    marker_yn::Vector{Int32} = model.markers.arrays.grid_marker_relationship.marker_yn.array
    marker_dx::Vector{Float64} = model.markers.arrays.grid_marker_relationship.marker_dx.array
    marker_dy::Vector{Float64} = model.markers.arrays.grid_marker_relationship.marker_dy.array

    marker_temperature::Vector{Float64} = model.markers.arrays.thermal.marker_TK.array

    tk1::Matrix{Float64} = model.heat_equation.arrays.temperature.tk1.array
    kt1::Matrix{Float64} = model.heat_equation.arrays.thermal_conductivity.kt1.array
    rhocp1::Matrix{Float64} = model.heat_equation.arrays.rhocp.rhocp1.array

    delta_temp_subgrid_markers = Vector{Float64}(undef, marknum)

    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == Int8(1)
            @inbounds begin
                # Define interpolation parameters
                dx_upr_left = marker_dx[imarker]
                dy_upr_left = marker_dy[imarker]
                ix_upr_left = marker_xn[imarker]
                iy_upr_left = marker_yn[imarker]
            end
            marker_temp_nodal = GridToMarker.get_marker_value_from_basic_grid_array(
                iy_upr_left, ix_upr_left, dy_upr_left, dx_upr_left, tk1)
            @inbounds delta_temp_subgrid_total = marker_temp_nodal - marker_temperature[imarker]
            marker_conductivity = GridToMarker.get_marker_value_from_basic_grid_array(
                iy_upr_left, ix_upr_left, dy_upr_left, dx_upr_left, kt1)
            marker_rhocp = GridToMarker.get_marker_value_from_basic_grid_array(
                iy_upr_left, ix_upr_left, dy_upr_left, dx_upr_left, rhocp1)
            diffusion_timescale = thermal_diffusion_timescale(
                marker_conductivity, marker_rhocp, xstpavg, ystpavg)
            subgrid_diffusion_term = calculate_subgrid_diffusion_term(
                timestep, diffusion_coef, diffusion_timescale)
            delta_temp_subgrid = delta_temp_subgrid_total * (1.0 - exp(subgrid_diffusion_term))
            @inbounds begin
                delta_temp_subgrid_markers[imarker] = delta_temp_subgrid
                # Correcting old temperature for the marker
                marker_temperature[imarker] += delta_temp_subgrid
            end
        else
            @inbounds delta_temp_subgrid_markers[imarker] = 0.0
        end
    end
    return delta_temp_subgrid_markers
end

function interpolate_subgrid_temperature_change_to_grid!(
    model::ModelData,
    delta_temp_subgrid_markers::Vector{Float64},
    inside_flags::Vector{Int8}
)::Nothing
    marknum::Int64 = model.markers.parameters.distribution.marknum.value

    marker_xn::Vector{Int32} = model.markers.arrays.grid_marker_relationship.marker_xn.array
    marker_yn::Vector{Int32} = model.markers.arrays.grid_marker_relationship.marker_yn.array

    marker_wtforULnode::Vector{Float64} = model.interpolation.arrays.marker_weights.marker_wtforULnode.array
    marker_wtforLLnode::Vector{Float64} = model.interpolation.arrays.marker_weights.marker_wtforLLnode.array
    marker_wtforURnode::Vector{Float64} = model.interpolation.arrays.marker_weights.marker_wtforURnode.array
    marker_wtforLRnode::Vector{Float64} = model.interpolation.arrays.marker_weights.marker_wtforLRnode.array

    fill!(model.heat_equation.arrays.temperature.dtkn.array, 0.0)
    dtkn::Matrix{Float64} = model.heat_equation.arrays.temperature.dtkn.array

    fill!(model.interpolation.arrays.grid_weights.wtnodes.array, 0.0)
    wtnodes::Matrix{Float64} = model.interpolation.arrays.grid_weights.wtnodes.array

    for imarker in 1:marknum
        if inside_flags[imarker] == Int8(1)
            @inbounds begin
                # Define interpolation parameters
                ix_upr_left = marker_xn[imarker]
                iy_upr_left = marker_yn[imarker]
                marker_weight_upr_left = marker_wtforULnode[imarker]
                marker_weight_lower_left = marker_wtforLLnode[imarker]
                marker_weight_upr_right = marker_wtforURnode[imarker]
                marker_weight_lower_right = marker_wtforLRnode[imarker]
                delta_temp_subgrid = delta_temp_subgrid_markers[imarker]
            end
            WeightFuncs.update_weight_inclusive!(
                iy_upr_left, ix_upr_left,
                marker_weight_upr_left, marker_weight_lower_left,
                marker_weight_upr_right, marker_weight_lower_right, 1.0,
                dtkn, delta_temp_subgrid
            )
            WeightFuncs.update_weight_inclusive!(
                iy_upr_left, ix_upr_left,
                marker_weight_upr_left, marker_weight_lower_left,
                marker_weight_upr_right, marker_weight_lower_right, 1.0,
                wtnodes, 1.0
            )
        end
    end
    BilinearAverage.calc_marker_average_subgrid_heat_diffusion!(model)
    return nothing
end

""" Remove subgrid effects from grid temperature changes.

This function updates the remaining grid conductive temperature change array `dtk1` 
by subtracting the subgrid temperature change array `dtkn` interpolated from 
markers.

# Arguments
- `model::ModelData`: The model data container

# Updated Arrays
## model.heat_equation.arrays.temperature
- dtk1.array: Array{Float64,2}
    - Remaining grid temperature change with subgrid effects removed.
"""
@inline function remove_subgrid_effects_from_grid_changes!(model::ModelData)
    dtk1 = model.heat_equation.arrays.temperature.dtk1
    dtkn = model.heat_equation.arrays.temperature.dtkn
    dtk1.array .= dtk1.array .- dtkn.array
end

""" Calculate the thermal diffusion timescale for subgrid diffusion.

# Arguments
- `marker_conductivity::Float64`: Thermal conductivity of the marker
- `marker_rhocp::Float64`: Density times heat capacity of the marker
- `xstpavg::Float64`: Average x-spacing of the basic grid
- `ystpavg::Float64`: Average y-spacing of the basic grid

# Returns
- `diffusion_timescale::Float64`: The thermal diffusion timescale for the marker
"""
@inline function thermal_diffusion_timescale(
    marker_conductivity::Float64,
    marker_rhocp::Float64,
    xstpavg::Float64,
    ystpavg::Float64
)::Float64
    diffusion_timescale = marker_rhocp / marker_conductivity / 
        (2.0 / xstpavg^2.0 + 2.0 / ystpavg^2.0)
    return diffusion_timescale
end

""" Calculate subgrid diffusion term.

The subgrid diffusion term appears in the exponential functions of the subgrid 
temperature change equation:

    dT_sg = dT_sg_total (1 - exp{-D_sg dt/t_sg)                        eq. 1

where dT_sg is the subgrid temperature change, dT_sg_total is the total subgrid 
temperature change, D_sg is the subgrid diffusion coefficient, dt is the 
timestep, and t_sg is the thermal diffusion timescale.

# Arguments
- `timestep::Float64`: Current timestep
- `diffusion_coef::Float64`: Subgrid diffusion coefficient
- `diffusion_timescale::Float64`: Thermal diffusion timescale

# Returns
- `subgrid_diffusion_term::Float64`: The subgrid diffusion term for the marker
"""
@inline function calculate_subgrid_diffusion_term(
    timestep::Float64,
    diffusion_coef::Float64,
    diffusion_timescale::Float64
)::Float64
    subgrid_diffusion_term = -diffusion_coef * timestep / diffusion_timescale
    subgrid_diffusion_term = max(subgrid_diffusion_term, -30.0)
    return subgrid_diffusion_term
end

end # module 