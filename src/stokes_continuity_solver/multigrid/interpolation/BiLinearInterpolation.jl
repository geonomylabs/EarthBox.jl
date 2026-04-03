module BiLinearInterpolation

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
    omx = 1.0 - dxUL
    omy = 1.0 - dyUL
    w00 = omx * omy
    w01 = omx * dyUL
    w10 = dxUL * omy
    w11 = dxUL * dyUL
    scalar_coarse[iULy  , jULx  ] += w00 * val
    weight_coarse[iULy  , jULx  ] += w00
    scalar_coarse[iULy+1, jULx  ] += w01 * val
    weight_coarse[iULy+1, jULx  ] += w01
    scalar_coarse[iULy  , jULx+1] += w10 * val
    weight_coarse[iULy  , jULx+1] += w10
    scalar_coarse[iULy+1, jULx+1] += w11 * val
    weight_coarse[iULy+1, jULx+1] += w11
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
    omx = 1.0 - dxUL
    omy = 1.0 - dyUL
    w00 = omx * omy
    w01 = omx * dyUL
    w10 = dxUL * omy
    w11 = dxUL * dyUL
    scalar_fine[i,j] += w00 * scalar_coarse[iULy  , jULx  ]
    scalar_fine[i,j] += w01 * scalar_coarse[iULy+1, jULx  ]
    scalar_fine[i,j] += w10 * scalar_coarse[iULy  , jULx+1]
    scalar_fine[i,j] += w11 * scalar_coarse[iULy+1, jULx+1]
    end # @inbounds
    return nothing
end

end