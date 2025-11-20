module WeightFuncs

"""
Functions used for calculating marker weight sums for bilinear interpolation.

The weight sums calculated in this module are used in the numerator of the
bilinear interpolation equation:

    grid_value_ij = SUM_m(marker_value_m*marker_weight_m) /
                    SUM_m(marker_weight_m)
"""

"""
    cell_volume_factor(xstp_b::Float64, ystp_b::Float64)::Float64

Make bilinear marker weight dependent on cell volume.

Marker weight is made inversely proportional to cell volume so if a marker
is in a relatively big cell it will have less of an impact.

# Arguments
- `xstp_b::Float64`: Width of grid cell in meters
- `ystp_b::Float64`: Height of grid cell in meters

# Returns
- `mwt_vol_fac::Float64`: Volume factor for marker weight
"""
function cell_volume_factor(xstp_b::Float64, ystp_b::Float64)::Float64
    ioption_local = 0
    if ioption_local == 0
        mwt_vol_fac = 1.0
    else
        mwt_vol_fac = 1.0/(xstp_b*ystp_b)
    end
    return mwt_vol_fac
end

""" Calculate weight sums for each node of grid cell that contains marker.

The weight sums are used in the numerator and denominator of the bilinear
interpolation equation:

    grid_value_ij = SUM_m(marker_value_m*marker_weight_m) /
                    SUM_m(marker_weight_m)

Set the marker scalar equal to 1.0 when calculating the weight sum for the
denominator.

This function includes any marker that is located a distance (dx_upr_left,
dy_upr_left) from the upper left node of the grid cell with indices
(iy_upr_left, ix_upr_left).

# Arguments
- `iy_upr_left::Int32`: Y-index of upper left node of grid cell
- `ix_upr_left::Int32`: X-index of upper left node of grid cell
- `marker_weight_upr_left::Float64`: Weight of marker at upper left node
- `marker_weight_lower_left::Float64`: Weight of marker at lower left node
- `marker_weight_upr_right::Float64`: Weight of marker at upper right node
- `marker_weight_lower_right::Float64`: Weight of marker at lower right node
- `vol_fac::Float64`: Volume factor for marker weight
- `weights::Matrix{Float64}`: Matrix of weights
- `marker_scalar::Float64`: Marker scalar

# Returns
- `nothing`: Nothing
"""
@inline function update_weight_inclusive!(
    iy_upr_left::Int32,
    ix_upr_left::Int32,
    marker_weight_upr_left::Float64,
    marker_weight_lower_left::Float64,
    marker_weight_upr_right::Float64,
    marker_weight_lower_right::Float64,
    vol_fac::Float64,
    weights::Matrix{Float64},
    marker_scalar::Float64
)::Nothing
    @inbounds begin
        # Upper left grid node
        weights[iy_upr_left, ix_upr_left] += marker_weight_upr_left*marker_scalar*vol_fac
        # Lower left grid node
        weights[iy_upr_left+1, ix_upr_left] += marker_weight_lower_left*marker_scalar*vol_fac
        # Upper right grid node
        weights[iy_upr_left, ix_upr_left+1] += marker_weight_upr_right*marker_scalar*vol_fac
        # Lower right grid node
        weights[iy_upr_left+1, ix_upr_left+1] += marker_weight_lower_right*marker_scalar*vol_fac
    end
    return nothing
end

""" Calculate weight sums for each node of grid cell that contains marker.

The weight sums are used in the numerator and denominator of the bilinear
interpolation equation:

    grid_value_ij = SUM_m(marker_value_m*marker_weight_m) /
                    SUM_m(marker_weight_m)

Set the marker scalar equal to 1.0 when calculating the weight sum for the
denominator.

This function includes any marker that is located a distance
(0.5dx_upr_left, 0.5dy_upr_left) from the upper left node of the grid cell
with indices (iy_upr_left, ix_upr_left).

# Arguments
- `iy_upr_left::Int32`: Y-index of upper left node of grid cell
- `ix_upr_left::Int32`: X-index of upper left node of grid cell
- `marker_weight_upr_left::Float64`: Weight of marker at upper left node
- `marker_weight_lower_left::Float64`: Weight of marker at lower left node
- `marker_weight_upr_right::Float64`: Weight of marker at upper right node
- `marker_weight_lower_right::Float64`: Weight of marker at lower right node
- `vol_fac::Float64`: Volume factor for marker weight
- `dx_upr_left::Float64`: X-direction distance to upper-left node of vy grid cell
- `dy_upr_left::Float64`: Y-direction distance to upper-left node of vy grid cell
- `weights::Matrix{Float64}`: Matrix of weights
- `marker_scalar::Float64`: Marker scalar

# Returns
- `nothing`: Nothing
"""
@inline function update_weight_exclusive!(
    iy_upr_left::Int32,
    ix_upr_left::Int32,
    marker_weight_upr_left::Float64,
    marker_weight_lower_left::Float64,
    marker_weight_upr_right::Float64,
    marker_weight_lower_right::Float64,
    vol_fac::Float64,
    dx_upr_left::Float64,
    dy_upr_left::Float64,
    weights::Matrix{Float64},
    marker_scalar::Float64
)::Nothing
    @inbounds begin
        # Upper left grid node
        if dx_upr_left <= 0.5 && dy_upr_left <= 0.5
            weights[iy_upr_left, ix_upr_left] += marker_weight_upr_left*marker_scalar*vol_fac
        end
        # Lower left grid node
        if dx_upr_left <= 0.5 <= dy_upr_left
            weights[iy_upr_left+1, ix_upr_left] += marker_weight_lower_left*marker_scalar*vol_fac
        end
        # Upper right grid node
        if dx_upr_left >= 0.5 >= dy_upr_left
            weights[iy_upr_left, ix_upr_left+1] += marker_weight_upr_right*marker_scalar*vol_fac
        end
        # Lower right grid node
        if dx_upr_left >= 0.5 && dy_upr_left >= 0.5
            weights[iy_upr_left+1, ix_upr_left+1] += marker_weight_lower_right*marker_scalar*vol_fac
        end
    end
    return nothing
end

""" Calculate weight sums for central node of grid cell that contains marker.

The weight sums are used in the numerator and denominator of the bilinear
interpolation equation:

    grid_value_ij = SUM_m(marker_value_m*marker_weight_m) /
                    SUM_m(marker_weight_m)

Set the marker scalar equal to 1.0 when calculating the weight sum for the
denominator.

# Arguments
- `iy_upr_left::Int32`: Y-index of upper left node of grid cell
- `ix_upr_left::Int32`: X-index of upper left node of grid cell
- `marker_weight_central::Float64`: Weight of marker at central node
- `vol_fac::Float64`: Volume factor for marker weight
- `weights::Matrix{Float64}`: Matrix of weights
- `marker_scalar::Float64`: Marker scalar

# Returns
- `nothing`: Nothing
"""
@inline function update_weight_central!(
    iy_upr_left::Int32,
    ix_upr_left::Int32,
    marker_weight_central::Float64,
    vol_fac::Float64,
    weights::Matrix{Float64},
    marker_scalar::Float64
)::Nothing
    @inbounds begin
        weights[iy_upr_left, ix_upr_left] += marker_weight_central*marker_scalar*vol_fac
    end
    return nothing
end

end # module 