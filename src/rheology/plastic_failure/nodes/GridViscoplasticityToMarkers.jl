module GridViscoplasticityToMarkers

using EarthBox.ModelDataContainer: ModelData
using EarthBox.GridFuncs: check_in_domain
using EarthBox.Interpolation: GridToMarker # To Do: Update once module is ready

function calc_marker_viscoplastic_viscosity!(
    model::ModelData, 
    inside_flags::Vector{Int8}
)::Nothing
    marker_xn = model.markers.arrays.grid_marker_relationship.marker_xn.array
    marker_yn = model.markers.arrays.grid_marker_relationship.marker_yn.array
    marker_dx = model.markers.arrays.grid_marker_relationship.marker_dx.array
    marker_dy = model.markers.arrays.grid_marker_relationship.marker_dy.array
    marker_pfailure = model.markers.arrays.rheology.marker_pfailure.array

    marknum = model.markers.parameters.distribution.marknum.value

    marker_eta = model.markers.arrays.rheology.marker_eta.array
    marker_eta_flow = model.markers.arrays.rheology.marker_eta_flow.array

    plastic_yield = model.stokes_continuity.arrays.plastic_def.plastic_yield.array
    etas1 = model.stokes_continuity.arrays.viscosity.etas1.array

    for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                xn_ul = marker_xn[imarker]
                yn_ul = marker_yn[imarker]
            end
            grid_yielding = check_basic_grid_for_yielding(plastic_yield, yn_ul, xn_ul)
            viscosity_flow = marker_eta_flow[imarker]
            if grid_yielding
                @inbounds begin
                    dx_upr_left = marker_dx[imarker]
                    dy_upr_left = marker_dy[imarker]
                end
                (
                    wt_upr_left_node, wt_lwr_left_node,
                    wt_upr_right_node, wt_lwr_right_node
                ) = GridToMarker.calc_marker_weights(dx_upr_left, dy_upr_left)
                viscosity_viscoplastic = interpolate_grid_viscosity_to_marker(
                    yn_ul, xn_ul,
                    plastic_yield, etas1,
                    wt_upr_left_node, wt_lwr_left_node,
                    wt_upr_right_node, wt_lwr_right_node
                    )
                if viscosity_viscoplastic >= viscosity_flow
                    @inbounds marker_eta[imarker] = viscosity_flow
                else
                    @inbounds begin
                        marker_eta[imarker] = viscosity_viscoplastic
                        marker_pfailure[imarker] = 1.0f0
                    end
                end
            else
                @inbounds marker_eta[imarker] = viscosity_flow
            end
        end
    end
    return nothing
end

@inline function check_basic_grid_for_yielding(
        plastic_yield::Array{Float64,2},
        iy_upr_left::Int32,
        ix_upr_left::Int32
)::Bool
    check = plastic_yield[iy_upr_left, ix_upr_left] > 0.0 ||
        plastic_yield[iy_upr_left+1, ix_upr_left] > 0.0 ||
        plastic_yield[iy_upr_left, ix_upr_left+1] > 0.0 ||
        plastic_yield[iy_upr_left+1, ix_upr_left+1] > 0.0
    return check
end

@inline function interpolate_grid_viscosity_to_marker(
    yn_ul::Int32,
    xn_ul::Int32,
    plastic_yield::Matrix{Float64},
    eta_bgrid::Matrix{Float64},
    wt_upr_left_node::Float64,
    wt_lwr_left_node::Float64,
    wt_upr_right_node::Float64,
    wt_lwr_right_node::Float64
)::Float64
    eta_ul = eta_bgrid[yn_ul, xn_ul]
    eta_ll = eta_bgrid[yn_ul+1, xn_ul]
    eta_ur = eta_bgrid[yn_ul, xn_ul+1]
    eta_lr = eta_bgrid[yn_ul+1, xn_ul+1]

    plastic_ul = plastic_yield[yn_ul, xn_ul]
    plastic_ll = plastic_yield[yn_ul+1, xn_ul]
    plastic_ur = plastic_yield[yn_ul, xn_ul+1]
    plastic_lr = plastic_yield[yn_ul+1, xn_ul+1]

    denom = (
        plastic_ul/eta_ul*wt_upr_left_node +
        plastic_ll/eta_ll*wt_lwr_left_node +
        plastic_ur/eta_ur*wt_upr_right_node +
        plastic_lr/eta_lr*wt_lwr_right_node
    )

    # Gerya 2019 used a numerator of 1 even though the sum of the weights is not 1 for the 4 nodes
    # if one of the plastic yield values is 0. This is not mathematically correct for the
    # harmonic average. To use the formulation of Gerya 2019 set use_sum_for_numerator to false.
    use_sum_for_numerator = true
    numerator = 1.0
    if use_sum_for_numerator
        numerator = (
            plastic_ul*wt_upr_left_node +
            plastic_ll*wt_lwr_left_node +
            plastic_ur*wt_upr_right_node +
            plastic_lr*wt_lwr_right_node
        )
    end

    return numerator/denom
end

end # module 