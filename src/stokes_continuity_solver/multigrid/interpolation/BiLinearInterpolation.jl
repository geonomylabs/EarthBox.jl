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
    @inbounds begin
    iULy  = grid_map.IULy[i,j]
    jULx  = grid_map.JULx[i,j]
    dyUL = grid_map.DyUL[i,j]
    dxUL = grid_map.DxUL[i,j] 
    val = scalar_fine[i,j]
    wt = WeightFuncs2d.upper_left_node_weight(dxUL, dyUL)
    scalar_coarse[iULy  , jULx  ] += wt * val
    weight_coarse[iULy  , jULx  ] += wt
    wt = WeightFuncs2d.lower_left_node_weight(dxUL, dyUL)
    scalar_coarse[iULy+1, jULx  ] += wt * val
    weight_coarse[iULy+1, jULx  ] += wt
    wt = WeightFuncs2d.upper_right_node_weight(dxUL, dyUL)
    scalar_coarse[iULy  , jULx+1] += wt * val
    weight_coarse[iULy  , jULx+1] += wt
    wt = WeightFuncs2d.lower_right_node_weight(dxUL, dyUL)
    scalar_coarse[iULy+1, jULx+1] += wt * val
    weight_coarse[iULy+1, jULx+1] += wt
    end # @inbounds
    return nothing
end

function coarse2fine_bilinear_interpolation!(
    i::Int,
    j::Int,
    grid_map::GridMapping2d,
    scalar_coarse::Array{Float64,2},
    scalar_fine::Array{Float64,2},
)::Nothing
    @inbounds begin
    iULy  = grid_map.IULy[i,j]
    jULx  = grid_map.JULx[i,j]
    dyUL = grid_map.DyUL[i,j]
    dxUL = grid_map.DxUL[i,j]
    wt = WeightFuncs2d.upper_left_node_weight(dxUL, dyUL)
    scalar_fine[i,j] += wt * scalar_coarse[iULy  , jULx  ]
    wt = WeightFuncs2d.lower_left_node_weight(dxUL, dyUL)
    scalar_fine[i,j] += wt * scalar_coarse[iULy+1, jULx  ]
    wt = WeightFuncs2d.upper_right_node_weight(dxUL, dyUL)
    scalar_fine[i,j] += wt * scalar_coarse[iULy  , jULx+1]
    wt = WeightFuncs2d.lower_right_node_weight(dxUL, dyUL)
    scalar_fine[i,j] += wt * scalar_coarse[iULy+1, jULx+1]
    end # @inbounds
    return nothing
end

end