module VelocityAndSpin

include("RungeKutta.jl")

import StaticArrays: MVector
import EarthBox.ModelDataContainer: ModelData
import EarthBox: GridFuncs
import EarthBox.Interpolation: MarkerGridMapping
import EarthBox.Interpolation: GridToMarker
import .RungeKutta

Base.@kwdef struct RKData
    marknum::Int64
    timestep::Float64
    xnum::Int
    ynum::Int
    marker_x::Vector{Float64}
    marker_y::Vector{Float64}
    gridx_b::Vector{Float64}
    gridy_b::Vector{Float64}
    gridx_vy::Vector{Float64}
    gridy_vx::Vector{Float64}
    xstp_b::Vector{Float64}
    ystp_b::Vector{Float64}
    xstp_vy::Vector{Float64}
    ystp_vx::Vector{Float64}
    marker_xn::Vector{Int32}
    marker_yn::Vector{Int32}
    vx1::Matrix{Float64}
    vy1::Matrix{Float64}
    esp::Matrix{Float64}
end

struct IterationData
    imarker::Int
    x_marker_initial::Float64
    y_marker_initial::Float64
end

""" Calculate velocity and spin for all markers using Runge-Kutta interpolation.

# Returns
- `Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}`: 
    (markers_vx, markers_vy, markers_spin)
"""
function interpolate_using_runge_kutta(
    model::ModelData,
    runge_kutta_order_max::Int,
    inside_flags::Vector{Int8}
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}

    rk_data = RKData(
        marknum = model.markers.parameters.distribution.marknum.value,
        timestep = model.timestep.parameters.main_time_loop.timestep.value,
        xnum = model.grids.parameters.geometry.xnum.value,
        ynum = model.grids.parameters.geometry.ynum.value,
        marker_x = model.markers.arrays.location.marker_x.array,
        marker_y = model.markers.arrays.location.marker_y.array,
        gridx_b = model.grids.arrays.basic.gridx_b.array,
        gridy_b = model.grids.arrays.basic.gridy_b.array,
        gridx_vy = model.grids.arrays.staggered_vy.gridx_vy.array,
        gridy_vx = model.grids.arrays.staggered_vx.gridy_vx.array,
        xstp_b = model.grids.arrays.basic.xstp_b.array,
        ystp_b = model.grids.arrays.basic.ystp_b.array,
        xstp_vy = model.grids.arrays.staggered_vy.xstp_vy.array,
        ystp_vx = model.grids.arrays.staggered_vx.ystp_vx.array,
        marker_xn = model.markers.arrays.grid_marker_relationship.marker_xn.array,
        marker_yn = model.markers.arrays.grid_marker_relationship.marker_yn.array,
        vx1 = model.stokes_continuity.arrays.staggered_grid_velocity.vx1.array,
        vy1 = model.stokes_continuity.arrays.staggered_grid_velocity.vy1.array,
        esp = model.stokes_continuity.arrays.strain_rate_and_spin.esp.array
    )

    markers_vx = Vector{Float64}(undef, rk_data.marknum)
    markers_vy = Vector{Float64}(undef, rk_data.marknum)
    markers_spin = Vector{Float64}(undef, rk_data.marknum)
    
    Threads.@threads for imarker in 1:rk_data.marknum
        if inside_flags[imarker] == Int8(1)
            @inbounds begin
                x_marker_initial = rk_data.marker_x[imarker]
                y_marker_initial = rk_data.marker_y[imarker]
            end
            iter_data = IterationData(imarker, x_marker_initial, y_marker_initial)

            x_marker_current = x_marker_initial
            y_marker_current = y_marker_initial

            vxm1 = vxm2 = vxm3 = vxm4 = 0.0
            vym1 = vym2 = vym3 = vym4 = 0.0
            espm1 = espm2 = espm3 = espm4 = 0.0
            for rk_order in 1:4
                vxm, vym, espm, x_marker_current, y_marker_current = rk_update( 
                    x_marker_current, y_marker_current, 
                    iter_data, rk_data, rk_order=rk_order
                )
                if rk_order == 1
                    vxm1, vym1, espm1 = vxm, vym, espm
                elseif rk_order == 2
                    vxm2, vym2, espm2 = vxm, vym, espm
                elseif rk_order == 3
                    vxm3, vym3, espm3 = vxm, vym, espm
                elseif rk_order == 4
                    vxm4, vym4, espm4 = vxm, vym, espm
                end
            end
            if runge_kutta_order_max == 4
                velocity_x = RungeKutta.apply_4th_order_runge_kutta(vxm1, vxm2, vxm3, vxm4)
                velocity_y = RungeKutta.apply_4th_order_runge_kutta(vym1, vym2, vym3, vym4)
                spin = RungeKutta.apply_4th_order_runge_kutta(espm1, espm2, espm3, espm4)
            else
                velocity_x = vxm1
                velocity_y = vym1
                spin = espm1
            end
            @inbounds begin
                markers_vx[imarker] = velocity_x
                markers_vy[imarker] = velocity_y
                markers_spin[imarker] = spin
            end
        else
            @inbounds begin
                markers_vx[imarker] = 0.0
                markers_vy[imarker] = 0.0
                markers_spin[imarker] = 0.0
            end
        end
    end
    return markers_vx, markers_vy, markers_spin
end

@inline function rk_update(
    x_marker_current::Float64,
    y_marker_current::Float64,
    iter_data::IterationData,
    rk_data::RKData;
    rk_order::Int64=1
)::Tuple{Float64, Float64, Float64, Float64, Float64}
    (
        x_index_upr_left_basic, y_index_upr_left_basic
    ) = get_upper_left_basic_grid_cell_indices(
        rk_order, iter_data.imarker, x_marker_current,
        y_marker_current, rk_data)
    (
        vxm, vym, espm
    ) = interpolate_for_velocity_and_spin(
        x_index_upr_left_basic, y_index_upr_left_basic,
        x_marker_current, y_marker_current, rk_data
    )
    x_marker_new, y_marker_new = advect_point(
        rk_order, iter_data.x_marker_initial, iter_data.y_marker_initial, 
        vxm, vym, rk_data.timestep
    )
    return vxm, vym, espm, x_marker_new, y_marker_new
end

@inline function get_upper_left_basic_grid_cell_indices(
    rk_order::Int,
    imarker::Int,
    x_marker_current::Float64,
    y_marker_current::Float64,
    rk_data::RKData,
)::Tuple{Int32, Int32}
    # If rk_order == 1, use the stored upper left grid cell indices for the marker
    if rk_order == 1
        @inbounds begin
            x_index_upr_left_basic = rk_data.marker_xn[imarker]
            y_index_upr_left_basic = rk_data.marker_yn[imarker]
        end
    else
        # Otherwise, use general grid mapping for an arbitrary x-y location
        x_index_upr_left_basic = MarkerGridMapping.get_index_of_left_node(
            x_marker_current, rk_data.gridx_b, rk_data.xnum)
        y_index_upr_left_basic = MarkerGridMapping.get_index_of_left_node(
            y_marker_current, rk_data.gridy_b, rk_data.ynum)
    end
    return x_index_upr_left_basic, y_index_upr_left_basic
end

@inline function interpolate_for_velocity_and_spin(
    x_index_upr_left_basic::Int32,
    y_index_upr_left_basic::Int32,
    x_marker_current::Float64,
    y_marker_current::Float64,
    rk_data::RKData,
)::Tuple{Float64, Float64, Float64}
    vxm = interpolate_vx(
        x_index_upr_left_basic, y_index_upr_left_basic,
        x_marker_current, y_marker_current,
        rk_data
    )
    vym = interpolate_vy(
        x_index_upr_left_basic, y_index_upr_left_basic,
        x_marker_current, y_marker_current,
        rk_data
    )
    espm = interpolate_spin(
        x_index_upr_left_basic, y_index_upr_left_basic,
        x_marker_current, y_marker_current,
        rk_data
    )
    return vxm, vym, espm
end

""" Calculate velocity for topography markers.

# Returns
- `Tuple{Vector{Float64}, Vector{Float64}}`: (markers_vx, markers_vy)
"""
function interpolate_velocity_for_topography_markers(
    model::ModelData,
    marker_x::Vector{Float64},
    marker_y::Vector{Float64},
    inside_flags::Vector{Int8},
    runge_kutta_order_max::Int
)::Tuple{Vector{Float64}, Vector{Float64}}
    marknum = length(marker_x)

    rk_data = RKData(
        marknum = marknum,
        timestep = model.timestep.parameters.main_time_loop.timestep.value,
        xnum = model.grids.parameters.geometry.xnum.value,
        ynum = model.grids.parameters.geometry.ynum.value,
        marker_x = marker_x,
        marker_y = marker_y,
        gridx_b = model.grids.arrays.basic.gridx_b.array,
        gridy_b = model.grids.arrays.basic.gridy_b.array,
        gridx_vy = model.grids.arrays.staggered_vy.gridx_vy.array,
        gridy_vx = model.grids.arrays.staggered_vx.gridy_vx.array,
        xstp_b = model.grids.arrays.basic.xstp_b.array,
        ystp_b = model.grids.arrays.basic.ystp_b.array,
        xstp_vy = model.grids.arrays.staggered_vy.xstp_vy.array,
        ystp_vx = model.grids.arrays.staggered_vx.ystp_vx.array,
        marker_xn = Vector{Int32}(undef, 1), # Not used
        marker_yn = Vector{Int32}(undef, 1), # Not used
        vx1 = model.stokes_continuity.arrays.staggered_grid_velocity.vx1.array,
        vy1 = model.stokes_continuity.arrays.staggered_grid_velocity.vy1.array,
        esp = Matrix{Float64}(undef, 1, 1) # Not used
    )

    markers_vx = zeros(Float64, marknum)
    markers_vy = zeros(Float64, marknum)

    for imarker in 1:marknum
        if inside_flags[imarker] == Int8(1)
            @inbounds begin
                x_marker_initial = marker_x[imarker]
                y_marker_initial = marker_y[imarker]
            end
            iter_data = IterationData(imarker, x_marker_initial, y_marker_initial)

            x_marker_current = x_marker_initial
            y_marker_current = y_marker_initial

            vxm1 = vxm2 = vxm3 = vxm4 = 0.0
            vym1 = vym2 = vym3 = vym4 = 0.0
            for rk_order in 1:4
                vxm, vym, x_marker_current, y_marker_current = rk_update_topography( 
                    x_marker_current, y_marker_current, 
                    iter_data, rk_data, rk_order=rk_order
                )
                if rk_order == 1
                    vxm1, vym1 = vxm, vym
                elseif rk_order == 2
                    vxm2, vym2 = vxm, vym
                elseif rk_order == 3
                    vxm3, vym3 = vxm, vym
                elseif rk_order == 4
                    vxm4, vym4 = vxm, vym
                end
            end
            if runge_kutta_order_max == 4
                velocity_x = RungeKutta.apply_4th_order_runge_kutta(vxm1, vxm2, vxm3, vxm4)
                velocity_y = RungeKutta.apply_4th_order_runge_kutta(vym1, vym2, vym3, vym4)
            else
                velocity_x = vxm1
                velocity_y = vym1
            end
            @inbounds begin
                markers_vx[imarker] = velocity_x
                markers_vy[imarker] = velocity_y
            end
        end
    end
    return markers_vx, markers_vy
end

@inline function rk_update_topography(
    x_marker_current::Float64,
    y_marker_current::Float64,
    iter_data::IterationData,
    rk_data::RKData;
    rk_order::Int64=1
)::Tuple{Float64, Float64, Float64, Float64}
    (
        x_index_upr_left_basic, y_index_upr_left_basic
    ) = get_upper_left_basic_grid_cell_indices_topography(
        x_marker_current, y_marker_current, rk_data)
    (
        vxm, vym
    ) = interpolate_for_velocity(
        x_index_upr_left_basic, y_index_upr_left_basic,
        x_marker_current, y_marker_current, rk_data
    )
    x_marker_new, y_marker_new = advect_point(
        rk_order, iter_data.x_marker_initial, iter_data.y_marker_initial, 
        vxm, vym, rk_data.timestep
    )
    return vxm, vym, x_marker_new, y_marker_new
end

@inline function get_upper_left_basic_grid_cell_indices_topography(
    x_marker_current::Float64,
    y_marker_current::Float64,
    rk_data::RKData,
)::Tuple{Int32, Int32}
    x_index_upr_left_basic = MarkerGridMapping.get_index_of_left_node(
        x_marker_current, rk_data.gridx_b, rk_data.xnum)
    y_index_upr_left_basic = MarkerGridMapping.get_index_of_left_node(
        y_marker_current, rk_data.gridy_b, rk_data.ynum)
    return x_index_upr_left_basic, y_index_upr_left_basic
end

@inline function interpolate_for_velocity(
    x_index_upr_left_basic::Int32,
    y_index_upr_left_basic::Int32,
    x_marker_current::Float64,
    y_marker_current::Float64,
    rk_data::RKData,
)::Tuple{Float64, Float64}
    vxm = interpolate_vx(
        x_index_upr_left_basic, y_index_upr_left_basic,
        x_marker_current, y_marker_current,
        rk_data
    )
    vym = interpolate_vy(
        x_index_upr_left_basic, y_index_upr_left_basic,
        x_marker_current, y_marker_current,
        rk_data
    )
    return vxm, vym
end

@inline function interpolate_vx(
    x_index_upr_left_basic::Int32,
    y_index_upr_left_basic::Int32,
    x_marker_current::Float64,
    y_marker_current::Float64,
    rk_data::RKData,
)::Float64
    # Interpolate velocity x
    x_index_upr_left, dx_upr_left = MarkerGridMapping.upr_left_x_mapping_vx_grid(
        x_index_upr_left_basic,
        x_marker_current,
        rk_data.gridx_b,
        rk_data.xstp_b
    )
    y_index_upr_left, dy_upr_left = MarkerGridMapping.upr_left_y_mapping_vx_grid(
        y_index_upr_left_basic,
        y_marker_current,
        rk_data.ynum,
        rk_data.gridy_vx,
        rk_data.ystp_vx
    )
    return GridToMarker.get_marker_value(
        y_index_upr_left,
        x_index_upr_left,
        dy_upr_left,
        dx_upr_left,
        rk_data.vx1
    )
end

@inline function interpolate_vy(
    x_index_upr_left_basic::Int32,
    y_index_upr_left_basic::Int32,
    x_marker_current::Float64,
    y_marker_current::Float64,
    rk_data::RKData,
)::Float64
    x_index_upr_left, dx_upr_left = MarkerGridMapping.upr_left_x_mapping_vy_grid(
        x_index_upr_left_basic,
        x_marker_current,
        rk_data.xnum,
        rk_data.gridx_vy,
        rk_data.xstp_vy
    )
    y_index_upr_left, dy_upr_left = MarkerGridMapping.upr_left_y_mapping_vy_grid(
        y_index_upr_left_basic,
        y_marker_current,
        rk_data.gridy_b,
        rk_data.ystp_b
    )
    return GridToMarker.get_marker_value(
        y_index_upr_left,
        x_index_upr_left,
        dy_upr_left,
        dx_upr_left,
        rk_data.vy1
    )
end

@inline function interpolate_spin(
    x_index_upr_left_basic::Int32,
    y_index_upr_left_basic::Int32,
    x_marker_current::Float64,
    y_marker_current::Float64,
    rk_data::RKData,
)::Float64
    x_index_upr_left, dx_upr_left = MarkerGridMapping.upr_left_x_mapping_basic_grid(
        x_index_upr_left_basic,
        x_marker_current,
        rk_data.gridx_b,
        rk_data.xstp_b
    )
    y_index_upr_left, dy_upr_left = MarkerGridMapping.upr_left_y_mapping_basic_grid(
        y_index_upr_left_basic,
        y_marker_current,
        rk_data.gridy_b,
        rk_data.ystp_b
    )
    return GridToMarker.get_marker_value(
        y_index_upr_left,
        x_index_upr_left,
        dy_upr_left,
        dx_upr_left,
        rk_data.esp
    )
end

@inline function advect_point(
    rk_order::Int,
    x_marker_initial::Float64,
    y_marker_initial::Float64,
    vxm::Float64,
    vym::Float64,
    timestep::Float64,
)::Tuple{Float64, Float64}
    # Advect point in velocity field
    x_marker_current = RungeKutta.calculate_x_for_next_runge_kutta_cycle(
        rk_order, timestep, x_marker_initial, vxm)
    y_marker_current = RungeKutta.calculate_y_for_next_runge_kutta_cycle(
        rk_order, timestep, y_marker_initial, vym)
    return x_marker_current, y_marker_current
end

function interpolate_velocity_for_topography_markers_old(
    model::ModelData,
    marker_x::Vector{Float64},
    marker_y::Vector{Float64},
    runge_kutta_order_max::Int
)::Tuple{Vector{Float64}, Vector{Float64}}
    marknum = length(marker_x)
    markers_vx = zeros(Float64, marknum)
    markers_vy = zeros(Float64, marknum)

    timestep = model.timestep.parameters.main_time_loop.timestep.value
    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value
    gridx_b = model.grids.arrays.basic.gridx_b.array
    gridy_b = model.grids.arrays.basic.gridy_b.array
    xstp_b = model.grids.arrays.basic.xstp_b.array
    ystp_b = model.grids.arrays.basic.ystp_b.array
    gridx_vy = model.grids.arrays.staggered_vy.gridx_vy.array
    gridy_vx = model.grids.arrays.staggered_vx.gridy_vx.array
    xstp_vy = model.grids.arrays.staggered_vy.xstp_vy.array
    ystp_vx = model.grids.arrays.staggered_vx.ystp_vx.array
    vx1 = model.stokes_continuity.arrays.staggered_grid_velocity.vx1.array
    vy1 = model.stokes_continuity.arrays.staggered_grid_velocity.vy1.array

    geometry = model.grids.parameters.geometry

    for imarker in 1:marknum
        vxm_rkarray = zeros(Float64, 4)
        vym_rkarray = zeros(Float64, 4)
        x_marker = marker_x[imarker]
        y_marker = marker_y[imarker]

        if GridFuncs.check_in_domain(geometry, x_marker, y_marker)
            x_marker_initial = marker_x[imarker]
            y_marker_initial = marker_y[imarker]
            x_marker_current = x_marker_initial
            y_marker_current = y_marker_initial

            for rk_order in 1:(runge_kutta_order_max)
                # Get indices of upper-left node
                x_index_upr_left_basic = MarkerGridMapping.get_index_of_left_node(
                    x_marker_current, gridx_b, xnum)
                y_index_upr_left_basic = MarkerGridMapping.get_index_of_left_node(
                    y_marker_current, gridy_b, ynum)

                # Interpolate velocity x
                x_index_upr_left, dx_upr_left = MarkerGridMapping.upr_left_x_mapping_vx_grid(
                    x_index_upr_left_basic,
                    x_marker_current,
                    gridx_b,
                    xstp_b
                )
                y_index_upr_left, dy_upr_left = MarkerGridMapping.upr_left_y_mapping_vx_grid(
                    y_index_upr_left_basic,
                    y_marker_current,
                    ynum,
                    gridy_vx,
                    ystp_vx
                )
                vxm_rkarray[rk_order] = GridToMarker.get_marker_value(
                    y_index_upr_left,
                    x_index_upr_left,
                    dy_upr_left,
                    dx_upr_left,
                    vx1
                )

                # Interpolate velocity y
                x_index_upr_left, dx_upr_left = MarkerGridMapping.upr_left_x_mapping_vy_grid(
                    x_index_upr_left_basic,
                    x_marker_current,
                    xnum,
                    gridx_vy,
                    xstp_vy
                )
                y_index_upr_left, dy_upr_left = MarkerGridMapping.upr_left_y_mapping_vy_grid(
                    y_index_upr_left_basic,
                    y_marker_current,
                    gridy_b,
                    ystp_b
                )
                vym_rkarray[rk_order] = GridToMarker.get_marker_value(
                    y_index_upr_left,
                    x_index_upr_left,
                    dy_upr_left,
                    dx_upr_left,
                    vy1
                )

                # Advect point in velocity field
                x_marker_current = RungeKutta.calculate_x_for_next_runge_kutta_cycle(
                    rk_order,
                    timestep,
                    x_marker_initial,
                    vxm_rkarray[rk_order]
                )
                y_marker_current = RungeKutta.calculate_y_for_next_runge_kutta_cycle(
                    rk_order,
                    timestep,
                    y_marker_initial,
                    vym_rkarray[rk_order]
                )
            end
            if runge_kutta_order_max == 4
                velocity_x = RungeKutta.apply_4th_order_runge_kutta(vxm_rkarray)
                velocity_y = RungeKutta.apply_4th_order_runge_kutta(vym_rkarray)
            else
                velocity_x = vxm_rkarray[1]
                velocity_y = vym_rkarray[1]
            end
            markers_vx[imarker] = velocity_x
            markers_vy[imarker] = velocity_y
        end
    end
    return markers_vx, markers_vy
end

end # module 