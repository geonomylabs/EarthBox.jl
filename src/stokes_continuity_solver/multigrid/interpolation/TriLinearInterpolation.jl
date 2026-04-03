module TriLinearInterpolation

import ..WeightFuncs3d
import ..GridMappingManager: GridMapping

""" Add to trilinear interpolation numerator and denominator

# Arguments
- `i::Int`:
    - Index of the node on the finer grid
- `j::Int`:
    - Index of the node on the finer grid
- `k::Int`:
    - Index of the node on the finer grid
- `grid_mapping::GridMapping`:
    - Grid mapping from finer to coarser grid
- `scalar_coarse::Array{Float64,3}`:
    - Coarser grid viscosity array
- `weight_coarse::Array{Float64,3}`:
    - Coarser grid weight array
- `scalar_fine::Array{Float64,3}`:
    - Filer grid viscosity array
"""
function add_to_numerator_and_denominator!(
    i::Int,
    j::Int,
    k::Int,
    grid_map::GridMapping,
    scalar_coarse::Array{Float64,3},
    weight_coarse::Array{Float64,3},
    scalar_fine::Array{Float64,3}
)::Nothing
    @inbounds begin
    iULBy  = grid_map.IULBy[i,j,k]
    jULBx  = grid_map.JULBx[i,j,k]
    kULBz  = grid_map.KULBz[i,j,k]
    dyULB = grid_map.DyULB[i,j,k]
    dxULB = grid_map.DxULB[i,j,k]
    dzULB = grid_map.DzULB[i,j,k]
    val = scalar_fine[i,j,k]
    wt = WeightFuncs3d.upper_left_back_node_weight(dxULB, dyULB, dzULB)
    scalar_coarse[iULBy  , jULBx  , kULBz  ] += wt * val
    weight_coarse[iULBy  , jULBx  , kULBz  ] += wt
    wt = WeightFuncs3d.upper_left_front_node_weight(dxULB, dyULB, dzULB)
    scalar_coarse[iULBy  , jULBx  , kULBz+1] += wt * val
    weight_coarse[iULBy  , jULBx  , kULBz+1] += wt
    wt = WeightFuncs3d.lower_left_back_node_weight(dxULB, dyULB, dzULB)
    scalar_coarse[iULBy+1, jULBx  , kULBz  ] += wt * val
    weight_coarse[iULBy+1, jULBx  , kULBz  ] += wt
    wt = WeightFuncs3d.lower_left_front_node_weight(dxULB, dyULB, dzULB)
    scalar_coarse[iULBy+1, jULBx  , kULBz+1] += wt * val
    weight_coarse[iULBy+1, jULBx  , kULBz+1] += wt
    wt = WeightFuncs3d.upper_right_back_node_weight(dxULB, dyULB, dzULB)
    scalar_coarse[iULBy  , jULBx+1, kULBz  ] += wt * val
    weight_coarse[iULBy  , jULBx+1, kULBz  ] += wt
    wt = WeightFuncs3d.upper_right_front_node_weight(dxULB, dyULB, dzULB)
    scalar_coarse[iULBy  , jULBx+1, kULBz+1] += wt * val
    weight_coarse[iULBy  , jULBx+1, kULBz+1] += wt
    wt = WeightFuncs3d.lower_right_back_node_weight(dxULB, dyULB, dzULB)
    scalar_coarse[iULBy+1, jULBx+1, kULBz  ] += wt * val
    weight_coarse[iULBy+1, jULBx+1, kULBz  ] += wt
    wt = WeightFuncs3d.lower_right_front_node_weight(dxULB, dyULB, dzULB)
    scalar_coarse[iULBy+1, jULBx+1, kULBz+1] += wt * val
    weight_coarse[iULBy+1, jULBx+1, kULBz+1] += wt
    end # @inbounds
    return nothing
end

function coarse2fine_trilinear_interpolation!(
    i::Int,
    j::Int,
    k::Int,
    grid_map::GridMapping,
    scalar_coarse::Array{Float64,3},
    scalar_fine::Array{Float64,3},
)::Nothing
    @inbounds begin
    iULBy  = grid_map.IULBy[i,j,k]
    jULBx  = grid_map.JULBx[i,j,k]
    kULBz  = grid_map.KULBz[i,j,k]
    dyULB = grid_map.DyULB[i,j,k]
    dxULB = grid_map.DxULB[i,j,k]
    dzULB = grid_map.DzULB[i,j,k]
    wt = WeightFuncs3d.upper_left_back_node_weight(dxULB, dyULB, dzULB)
    scalar_fine[i,j,k] += wt * scalar_coarse[iULBy  , jULBx  , kULBz  ]
    wt = WeightFuncs3d.upper_left_front_node_weight(dxULB, dyULB, dzULB)
    scalar_fine[i,j,k] += wt * scalar_coarse[iULBy  , jULBx  , kULBz+1]
    wt = WeightFuncs3d.lower_left_back_node_weight(dxULB, dyULB, dzULB)
    scalar_fine[i,j,k] += wt * scalar_coarse[iULBy+1, jULBx  , kULBz  ]
    wt = WeightFuncs3d.lower_left_front_node_weight(dxULB, dyULB, dzULB)
    scalar_fine[i,j,k] += wt * scalar_coarse[iULBy+1, jULBx  , kULBz+1]
    wt = WeightFuncs3d.upper_right_back_node_weight(dxULB, dyULB, dzULB)
    scalar_fine[i,j,k] += wt * scalar_coarse[iULBy  , jULBx+1, kULBz  ]
    wt = WeightFuncs3d.upper_right_front_node_weight(dxULB, dyULB, dzULB)
    scalar_fine[i,j,k] += wt * scalar_coarse[iULBy  , jULBx+1, kULBz+1]
    wt = WeightFuncs3d.lower_right_back_node_weight(dxULB, dyULB, dzULB)
    scalar_fine[i,j,k] += wt * scalar_coarse[iULBy+1, jULBx+1, kULBz  ]
    wt = WeightFuncs3d.lower_right_front_node_weight(dxULB, dyULB, dzULB)
    scalar_fine[i,j,k] += wt * scalar_coarse[iULBy+1, jULBx+1, kULBz+1]
    end # @inbounds
    return nothing
end

end