module TriLinearInterpolation

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
    # Fused trilinear weights: w = (1-dx|dx)*(1-dy|dy)*(1-dz|dz) over the 8 corners.
    omx = 1.0 - dxULB
    omy = 1.0 - dyULB
    omz = 1.0 - dzULB
    w000 = omx * omy * omz
    w001 = omx * omy * dzULB
    w010 = omx * dyULB * omz
    w011 = omx * dyULB * dzULB
    w100 = dxULB * omy * omz
    w101 = dxULB * omy * dzULB
    w110 = dxULB * dyULB * omz
    w111 = dxULB * dyULB * dzULB
    scalar_coarse[iULBy  , jULBx  , kULBz  ] += w000 * val
    weight_coarse[iULBy  , jULBx  , kULBz  ] += w000
    scalar_coarse[iULBy  , jULBx  , kULBz+1] += w001 * val
    weight_coarse[iULBy  , jULBx  , kULBz+1] += w001
    scalar_coarse[iULBy+1, jULBx  , kULBz  ] += w010 * val
    weight_coarse[iULBy+1, jULBx  , kULBz  ] += w010
    scalar_coarse[iULBy+1, jULBx  , kULBz+1] += w011 * val
    weight_coarse[iULBy+1, jULBx  , kULBz+1] += w011
    scalar_coarse[iULBy  , jULBx+1, kULBz  ] += w100 * val
    weight_coarse[iULBy  , jULBx+1, kULBz  ] += w100
    scalar_coarse[iULBy  , jULBx+1, kULBz+1] += w101 * val
    weight_coarse[iULBy  , jULBx+1, kULBz+1] += w101
    scalar_coarse[iULBy+1, jULBx+1, kULBz  ] += w110 * val
    weight_coarse[iULBy+1, jULBx+1, kULBz  ] += w110
    scalar_coarse[iULBy+1, jULBx+1, kULBz+1] += w111 * val
    weight_coarse[iULBy+1, jULBx+1, kULBz+1] += w111
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
    omx = 1.0 - dxULB
    omy = 1.0 - dyULB
    omz = 1.0 - dzULB
    w000 = omx * omy * omz
    w001 = omx * omy * dzULB
    w010 = omx * dyULB * omz
    w011 = omx * dyULB * dzULB
    w100 = dxULB * omy * omz
    w101 = dxULB * omy * dzULB
    w110 = dxULB * dyULB * omz
    w111 = dxULB * dyULB * dzULB
    scalar_fine[i,j,k] += w000 * scalar_coarse[iULBy  , jULBx  , kULBz  ]
    scalar_fine[i,j,k] += w001 * scalar_coarse[iULBy  , jULBx  , kULBz+1]
    scalar_fine[i,j,k] += w010 * scalar_coarse[iULBy+1, jULBx  , kULBz  ]
    scalar_fine[i,j,k] += w011 * scalar_coarse[iULBy+1, jULBx  , kULBz+1]
    scalar_fine[i,j,k] += w100 * scalar_coarse[iULBy  , jULBx+1, kULBz  ]
    scalar_fine[i,j,k] += w101 * scalar_coarse[iULBy  , jULBx+1, kULBz+1]
    scalar_fine[i,j,k] += w110 * scalar_coarse[iULBy+1, jULBx+1, kULBz  ]
    scalar_fine[i,j,k] += w111 * scalar_coarse[iULBy+1, jULBx+1, kULBz+1]
    end # @inbounds
    return nothing
end

end