module BiLinearInterpolation

import ..WeightFuncs2d
import ..GridMappingManager: GridMapping2d

""" Add to bilinear interpolation numerator and denominator

# Arguments
- `i::Int`:
    - Index of the node on the finer grid
- `j::Int`:
    - Index of the node on the finer grid
- `grid_mapping::GridMapping2d`:
    - Grid mapping from finer to coarser grid
- `scalar_coarse::Array{Float64,2}`:
    - Coarser grid viscosity array
- `weight_coarse::Array{Float64,2}`:
    - Coarser grid weight array
- `scalar_fine::Array{Float64,2}`:
    - Filer grid viscosity array
"""
function add_to_numerator_and_denominator_2d!(
    i::Int,
    j::Int,
    grid_map::GridMapping2d,
    scalar_coarse::Array{Float64,2},
    weight_coarse::Array{Float64,2},
    scalar_fine::Array{Float64,2}
)::Nothing
    # Indices of upper-left node of coarse grid cell containing the finer grid node (i,j)
    iULy  = grid_map.IULy[i,j]
    jULx  = grid_map.JULx[i,j]
    # Normalized distance between finer grid node (i,j) and upper-left node of coarse grid cell
    dyUL = grid_map.DyUL[i,j]
    dxUL = grid_map.DxUL[i,j] 
    # Add scalar on finer grid to the coarser grid for 4 nodes bounding the cell
    # Upper-left reference node
    wt = WeightFuncs2d.upper_left_node_weight(dxUL, dyUL)
    scalar_coarse[iULy  , jULx  ] += wt * scalar_fine[i,j]
    weight_coarse[iULy  , jULx  ] += wt
    # Lower-left node
    wt = WeightFuncs2d.lower_left_node_weight(dxUL, dyUL)
    scalar_coarse[iULy+1, jULx  ] += wt * scalar_fine[i,j]
    weight_coarse[iULy+1, jULx  ] += wt
    # Upper-right node
    wt = WeightFuncs2d.upper_right_node_weight(dxUL, dyUL)
    scalar_coarse[iULy  , jULx+1] += wt * scalar_fine[i,j]
    weight_coarse[iULy  , jULx+1] += wt
    # Lower-right node
    wt = WeightFuncs2d.lower_right_node_weight(dxUL, dyUL)
    scalar_coarse[iULy+1, jULx+1] += wt * scalar_fine[i,j]
    weight_coarse[iULy+1, jULx+1] += wt
    return nothing
end

function coarse2fine_bilinear_interpolation!(
    i::Int,
    j::Int,
    grid_map::GridMapping2d,
    scalar_coarse::Array{Float64,2},
    scalar_fine::Array{Float64,2},
)::Nothing
    # Indices of upper-left node of coarse grid cell containing the finer grid node (i,j)
    iULy  = grid_map.IULy[i,j]
    jULx  = grid_map.JULx[i,j]
    # Normalized distance between finer grid node (i,j) and upper-left node of coarse grid cell
    dyUL = grid_map.DyUL[i,j]
    dxUL = grid_map.DxUL[i,j]
    # Add scalar from the coarser to finer level 
    # using bilinear interpolation from 4 nodes bounding the cell
    wt = WeightFuncs2d.upper_left_node_weight(dxUL, dyUL)
    scalar_fine[i,j] += wt * scalar_coarse[iULy  , jULx  ]
    wt = WeightFuncs2d.lower_left_node_weight(dxUL, dyUL)
    scalar_fine[i,j] += wt * scalar_coarse[iULy+1, jULx  ]
    wt = WeightFuncs2d.upper_right_node_weight(dxUL, dyUL)
    scalar_fine[i,j] += wt * scalar_coarse[iULy  , jULx+1]
    wt = WeightFuncs2d.lower_right_node_weight(dxUL, dyUL)
    scalar_fine[i,j] += wt * scalar_coarse[iULy+1, jULx+1]
    return nothing
end

end