"""
    MarkerGridNodeWeights

Module containing functions to calculate marker-to-grid interpolation parameters.

Marker information is interpolated to grid nodes using first-order bilinear
weight factors that are formulated using the normalized x and y-distances (i.e.
dx and dy) between the upper-left node and the marker. Marker weight-factor
formulas are as follows:

    mwtUL = (1.0 - dx)*(1.0 - dy)
    mwtLL = (1.0 - dx)*dy
    mwtUR = dx*(1.0 - dy)
    mwtLR = dx*dy
    mwtC = (1.0 - abs(0.5 - dx))*(1.0 - abs(0.5 - dy))

where mwt refers to marker weight and UL, LL, UR, LR and C refer to the
upper-left node, lower-left node, upper right-node, lower-right node and
central node respectively.

This upper-left node based formulation enables the relationship between nodes
and markers to be quickly computed using the pre-computed upper-left node index
for each marker.

                 j  xn                           xn+1

        i  yn    grid(yn,xn)--------------------grid(yn,xn+1)
                    |           ^                  |
                    |           |                  |
                    |          dy                  |
                    |           |                  |
                    |           v                  |
                    |<----dx--->x marker(ym,xm)    |
                    |                              |
                    |                              |
           yn+1  grid(yn+1,xn)-------------------grid(yn+1, xn+1)

    Figure 1: Schematic showing how surrounding grid nodes are related
    to a marker located within a grid cell. The upper-left grid node
    (yn, xn) is calculated for a given marker location (ym, xm).
"""
module MarkerGridNodeWeights

using EarthBox.ModelDataContainer: ModelData
using EarthBox: GridFuncs

"""
    calc_marker_weights_for_basic_and_pressure_grids!(model::ModelData)

Calculate weights for each marker for first-order bilinear scheme.

# Arguments
- `model`: Model data structure containing marker and grid information

# Updated Arrays
## Updated arrays from group `model.interpolation.arrays.marker_weights`
- `marker_wtforULnode.array::Vector{Float64}` (marknum):
    - Marker weight for upper-left node
- `marker_wtforLLnode.array::Vector{Float64}` (marknum):
    - Marker weight for lower-left node
- `marker_wtforURnode.array::Vector{Float64}` (marknum):
    - Marker weight for upper-right node
- `marker_wtforLRnode.array::Vector{Float64}` (marknum):
    - Marker weight for lower-right node
- `marker_wtforCnode.array::Vector{Float64}` (marknum):
    - Marker weight for central node
"""
function calc_marker_weights_for_basic_and_pressure_grids!(
    model::ModelData,
    inside_flags::Vector{Int8}
)::Nothing
    marknum = model.markers.parameters.distribution.marknum.value

    marker_dx = model.markers.arrays.grid_marker_relationship.marker_dx.array
    marker_dy = model.markers.arrays.grid_marker_relationship.marker_dy.array

    marker_weights = model.interpolation.arrays.marker_weights
    marker_wtforULnode = marker_weights.marker_wtforULnode.array
    marker_wtforLLnode = marker_weights.marker_wtforLLnode.array
    marker_wtforURnode = marker_weights.marker_wtforURnode.array
    marker_wtforLRnode = marker_weights.marker_wtforLRnode.array
    marker_wtforCnode = marker_weights.marker_wtforCnode.array

    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            dx_upper_left = marker_dx[imarker]
            dy_upper_left = marker_dy[imarker]
            # Upper-left node
            marker_wtforULnode[imarker] = (1.0 - dx_upper_left) * (1.0 - dy_upper_left)
            # Lower-left node
            marker_wtforLLnode[imarker] = (1.0 - dx_upper_left) * dy_upper_left
            # Upper-right node
            marker_wtforURnode[imarker] = dx_upper_left * (1.0 - dy_upper_left)
            # Lower-right node
            marker_wtforLRnode[imarker] = dx_upper_left * dy_upper_left
            # Central node
            marker_wtforCnode[imarker] = (1.0 - abs(0.5 - dx_upper_left)) * (1.0 - abs(0.5 - dy_upper_left))
        end
    end
    return nothing
end

"""
    calc_marker_weights_for_vy_grid!(model::ModelData)

Calculate weights for each marker for first-order bilinear scheme on Vy grid.

# Arguments
- `model`: Model data structure containing marker and grid information

# Updated Arrays
## Updated arrays from group `model.interpolation.arrays.marker_weights`
- `marker_wtforULnodeVy.array::Vector{Float64}` (marknum):
    - Marker weight for upper-left node of staggered Vy grid
- `marker_wtforLLnodeVy.array::Vector{Float64}` (marknum):
    - Marker weight for lower-left node of staggered Vy grid
- `marker_wtforURnodeVy.array::Vector{Float64}` (marknum):
    - Marker weight for upper-right node of staggered Vy grid
- `marker_wtforLRnodeVy.array::Vector{Float64}` (marknum):
    - Marker weight for lower-right node of staggered Vy grid
"""
function calc_marker_weights_for_vy_grid!(model::ModelData, inside_flags::Vector{Int8})
    marknum = model.markers.parameters.distribution.marknum.value

    marker_dx_vy = model.markers.arrays.grid_marker_relationship.marker_dx_vy.array
    marker_dy_vy = model.markers.arrays.grid_marker_relationship.marker_dy_vy.array

    marker_weights = model.interpolation.arrays.marker_weights
    marker_wtforULnodeVy = marker_weights.marker_wtforULnodeVy.array
    marker_wtforLLnodeVy = marker_weights.marker_wtforLLnodeVy.array
    marker_wtforURnodeVy = marker_weights.marker_wtforURnodeVy.array
    marker_wtforLRnodeVy = marker_weights.marker_wtforLRnodeVy.array

    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            dx_upper_left = marker_dx_vy[imarker]
            dy_upper_left = marker_dy_vy[imarker]
            
            # Upper-left node
            marker_wtforULnodeVy[imarker] = (1.0 - dx_upper_left) * (1.0 - dy_upper_left)
            # Lower-left node
            marker_wtforLLnodeVy[imarker] = (1.0 - dx_upper_left) * dy_upper_left
            # Upper-right node
            marker_wtforURnodeVy[imarker] = dx_upper_left * (1.0 - dy_upper_left)
            # Lower-right node
            marker_wtforLRnodeVy[imarker] = dx_upper_left * dy_upper_left
        end
    end
    return nothing
end

end # module 