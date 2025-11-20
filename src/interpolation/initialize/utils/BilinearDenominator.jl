module BilinearDenominator

import EarthBox.ModelDataContainer.ModelData
import EarthBox.ModelDataContainer.InterpolationContainer: clear_grid_weights!
import EarthBox.GridFuncs: check_in_domain
import ...WeightFuncs

"""
    calc_marker_weight_sums_for_basic_and_pressure_nodes(model::ModelData)

Calculate marker weight sums for denominator at grid nodes.

Marker weight factors are summed at grid nodes for the denominator of the
first-order bilinear interpolation equation.

# Updated Arrays
## Updated arrays from group model.interpolation.arrays.grid_weights
- `wtnodes::Matrix{Float64}` (ynum, xnum):
    - Summed marker weight factors at basic grid nodes used in
    denominator of average for interpolating marker density,
    temperature, conductivity, heat capacity, and heat production to
    basic grid nodes.
- `wtetas::Matrix{Float64}` (ynum, xnum):
    - Summed marker weight factors at basic grid nodes used in
    denominator of average for interpolating marker shear viscosity,
    shear modulus under shear stress, shear stress, and plastic flags
    to basic grid nodes.
- `wtetan::Matrix{Float64}` (ynum-1, xnum-1):
    - Summed marker weight factors at pressure grid nodes used in
    denominator of average for interpolating marker normal viscosity,
    shear modulus under normal stress, normal stress, and plastic flags
    to pressure grid nodes.
"""
function calc_marker_weight_sums_for_basic_and_pressure_nodes!(model::ModelData, inside_flags::Vector{Int8})
    # Number of markers
    marknum = model.markers.parameters.distribution.marknum.value
    # X- and Y-spacing of basic grid
    #xstp_b = model.grids.arrays.basic.xstp_b.array
    #ystp_b = model.grids.arrays.basic.ystp_b.array
    # X- and Y indices of upper left node of grid cell
    marker_xn = model.markers.arrays.grid_marker_relationship.marker_xn.array
    marker_yn = model.markers.arrays.grid_marker_relationship.marker_yn.array
    # X- and Y distances from upper left node of grid cell
    marker_dx = model.markers.arrays.grid_marker_relationship.marker_dx.array
    marker_dy = model.markers.arrays.grid_marker_relationship.marker_dy.array
    # Distance-based marker weights for each node of grid cell
    marker_wtforULnode = model.interpolation.arrays.marker_weights.marker_wtforULnode.array
    marker_wtforLLnode = model.interpolation.arrays.marker_weights.marker_wtforLLnode.array
    marker_wtforURnode = model.interpolation.arrays.marker_weights.marker_wtforURnode.array
    marker_wtforLRnode = model.interpolation.arrays.marker_weights.marker_wtforLRnode.array
    marker_wtforCnode = model.interpolation.arrays.marker_weights.marker_wtforCnode.array
    # Marker weight sums at grid nodes
    wtnodes = model.interpolation.arrays.grid_weights.wtnodes.array
    wtetas = model.interpolation.arrays.grid_weights.wtetas.array
    wtetan = model.interpolation.arrays.grid_weights.wtetan.array

    for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                # Upper left x- and y- indices of grid cell
                ix_upr_left = marker_xn[imarker]
                iy_upr_left = marker_yn[imarker]
                # X- and Y-distances from upper left node of grid cell
                dx_upr_left = marker_dx[imarker]
                dy_upr_left = marker_dy[imarker]
                # Distance-based marker weight factors for each node of grid cell
                marker_weight_upr_left = marker_wtforULnode[imarker]
                marker_weight_lower_left = marker_wtforLLnode[imarker]
                marker_weight_upr_right = marker_wtforURnode[imarker]
                marker_weight_lower_right = marker_wtforLRnode[imarker]
                marker_weight_central = marker_wtforCnode[imarker]
            end
            # Volume factor for marker weight factors
            #vol_fac = WeightFuncs.cell_volume_factor(xstp_b[ix_upr_left], ystp_b[iy_upr_left])
            # Update marker weight sum for each node of basic grid using the
            # inclusive method.
            WeightFuncs.update_weight_inclusive!(
                iy_upr_left, ix_upr_left,
                marker_weight_upr_left, marker_weight_lower_left,
                marker_weight_upr_right, marker_weight_lower_right, 1.0,
                wtnodes, 1.0
            )
            # Update marker weight sum for central node of basic grid using the
            # exclusive method
            WeightFuncs.update_weight_exclusive!(
                iy_upr_left, ix_upr_left,
                marker_weight_upr_left, marker_weight_lower_left,
                marker_weight_upr_right, marker_weight_lower_right, 1.0,
                dx_upr_left, dy_upr_left,
                wtetas, 1.0
            )
            # Update marker weight sum for each node of pressure grid using the
            # inclusive method
            WeightFuncs.update_weight_central!(
                iy_upr_left, ix_upr_left, marker_weight_central, 1.0,
                wtetan, 1.0
            )
        end
    end
end

"""
    calc_marker_weight_sums_for_vy_nodes(model::ModelData)

Calculate marker weight sums for denominator at vy grid nodes.

Marker weight factors are summed at grid nodes for the denominator of the
first-order bilinear interpolation equation.

# Updated Arrays
## Updated arrays from group model.interpolation.arrays.grid_weights
- `wtnodes_vy::Matrix{Float64}` (ynum, xnum+1):
    - Summed marker weight factors at Vy grid nodes used in denominator
    of average for interpolating density to Vy grid nodes.
"""
function calc_marker_weight_sums_for_vy_nodes!(model::ModelData, inside_flags::Vector{Int8})
    # Number of markers
    marknum = model.markers.parameters.distribution.marknum.value
    # X- and Y-spacing of vy grid
    #xstp_vy = model.grids.arrays.staggered_vy.xstp_vy.array
    #ystp_b = model.grids.arrays.basic.ystp_b.array
    # X- and Y indices of upper left node of vy grid cell
    marker_xn_vy = model.markers.arrays.grid_marker_relationship.marker_xn_vy.array
    marker_yn_vy = model.markers.arrays.grid_marker_relationship.marker_yn_vy.array
    # Distance-based marker weights for each node of grid cell
    marker_wtforULnodeVy = model.interpolation.arrays.marker_weights.marker_wtforULnodeVy.array
    marker_wtforLLnodeVy = model.interpolation.arrays.marker_weights.marker_wtforLLnodeVy.array
    marker_wtforURnodeVy = model.interpolation.arrays.marker_weights.marker_wtforURnodeVy.array
    marker_wtforLRnodeVy = model.interpolation.arrays.marker_weights.marker_wtforLRnodeVy.array
    # Marker weight sums at grid nodes
    wtnodes_vy = model.interpolation.arrays.grid_weights.wtnodes_vy.array

    for imarker in 1:marknum
        @inbounds begin
            if inside_flags[imarker] == 1
                # Upper left x- and y- indices of grid cell
                ix_upr_left = marker_xn_vy[imarker]
                iy_upr_left = marker_yn_vy[imarker]
                # Distance-based marker weight factors for each node of grid cell
                marker_weight_upr_left = marker_wtforULnodeVy[imarker]
                marker_weight_lower_left = marker_wtforLLnodeVy[imarker]
                marker_weight_upr_right = marker_wtforURnodeVy[imarker]
                marker_weight_lower_right = marker_wtforLRnodeVy[imarker]
                # Volume factor for marker weight factors
                #vol_fac = WeightFuncs.cell_volume_factor(xstp_vy[ix_upr_left], ystp_b[iy_upr_left])
                # Update marker weight sum for each node of basic grid using the
                # inclusive method.
                WeightFuncs.update_weight_inclusive!(
                    iy_upr_left, ix_upr_left,
                    marker_weight_upr_left, marker_weight_lower_left,
                    marker_weight_upr_right, marker_weight_lower_right, 1.0,
                    wtnodes_vy, 1.0
                )
            end
        end
    end
end

"""
    clear_marker_weight_sums_on_grids(model::ModelData)

Clear marker weight sums for denominator or bilinear interp. equation.

# Updated Numpy Arrays
## Updated arrays from group model.interpolation.arrays.grid_weights
- `wtnodes::Matrix{Float64}` (ynum, xnum):
    - Summed marker weight factors at basic grid nodes used in denominator of
    average for interpolating marker density, temperature, conductivity, 
    heat capacity, and heat production to basic grid nodes.
- `wtetas::Matrix{Float64}` (ynum, xnum):
    - Summed marker weight factors at basic grid nodes used in
    denominator of average for interpolating marker shear viscosity,
    shear modulus under shear stress, shear stress, and plastic flags
    to basic grid nodes.
- `wtetan::Matrix{Float64}` (ynum, xnum):
    - Summed marker weight factors at pressure grid nodes used in
    denominator of average for interpolating marker normal viscosity,
    shear modulus under normal stress, normal stress, and plastic flags
    to pressure grid nodes.
- `wtnodes_vy::Matrix{Float64}` (ynum, xnum):
    - Summed marker weight factors at Vy grid nodes used in denominator
    of average for interpolating density to Vy grid nodes.
"""
function clear_marker_weight_sums_on_grids!(model::ModelData)
    clear_grid_weights!(model.interpolation)
end

end # module 