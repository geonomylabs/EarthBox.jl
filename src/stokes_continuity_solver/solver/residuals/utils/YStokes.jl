"""
    YStokes

Compute residual for y-Stokes equation at a given grid location (i,j).

+-------------------- -+-------vy(i-1,j)------+----------------------+-----
|                      |                      |                      |
|                      |                      |                      |
|                  vx(i,j-1)  P(i-1,j-1)   vx(i,j)                   |ystp(i-1)-------
|                      |      etan(i-1,j-1)   |                      |
|                      |                      |                      |
+-----vy(i,j-1)---etas(i,j-1)----vy(i,j)----etas(i,j)---vy(i,j+1)----+----- ystpc(i)
|                      |                      |                      |
|                      |                      |                      |
|                  vx(i+1,j-1)  P(i,j-1)   vx(i+1,j)                 |ystp(i)------
|                      |      etan(i,j-1)     |                      |
|                      |                      |                      |
+----------------------+-------vy(i+1,j)------+----------------------+-----
          |        xstpc(j-1)      |        xstpc(j)      |
                        |       xstp(j-1)      |

Figure 1: Y-Stokes equation stencil used to calculate y-Stokes residual.
"""
module YStokes

using EarthBox.MathTools
using ..ResidualStructs: StokesResidualsInput

"""
    calculate_residual!(
        residuals_input::StokesResidualsInput,
        i::Int,
        j::Int,
        resy::Matrix{Float64}
    )::Nothing

Compute residual for y-Stokes equation at a given grid location (i,j).

The y-Stokes equation is defined as follows:

    dsyy/dy + dsyx/dx - dP/dy = RY

The residual after solving the system of equations is as follows:

    RY - dsyy/dy - dsyx/dx + dP/dy = residual

See Figure 1 above for details on the stencil used to calculate
derivatives.

# Arguments
- `residuals_input::StokesResidualsInput`: Input parameters for residual calculation
- `i::Int`: Grid index in y direction
- `j::Int`: Grid index in x direction
- `resy::Matrix{Float64}`: Residual array for y-Stokes equation

# Returns
- `Nothing`: Updates resy array in place
"""
function calculate_residual!(
    residuals_input::StokesResidualsInput,
    i::Int,
    j::Int,
    resy::Matrix{Float64}
)::Nothing
    xnum = residuals_input.xnum
    ynum = residuals_input.ynum
    xstp_b = residuals_input.xstp_b
    ystp_b = residuals_input.ystp_b
    xstp_vy = residuals_input.xstp_vy
    ystp_vx = residuals_input.ystp_vx
    bintern_zone = residuals_input.bintern_zone
    etas0 = residuals_input.etas0
    etan0 = residuals_input.etan0
    RY1 = residuals_input.RY1
    vx1 = residuals_input.vx1
    vy1 = residuals_input.vy1
    pr1 = residuals_input.pr1

    if is_valid_node_for_y_stokes(i, j, ynum, bintern_zone)
        # vy-boundary conditions
        if is_boundary_node(i, j, xnum, ynum)
            resy[i,j] = 0.0
        else
            # dsyy/dy-dP/dy
            resy[i,j] = (
                RY1[i,j]
                - (
                    2.0 * (
                        etan0[i,j-1] * (vy1[i+1,j] - vy1[i,j]) / ystp_b[i]
                        - etan0[i-1,j-1] * (vy1[i,j] - vy1[i-1,j]) / ystp_b[i-1]
                    )
                    - (pr1[i,j-1] - pr1[i-1,j-1])
                ) / ystp_vx[i]
            )
            # dsxy/dx
            resy[i,j] = (
                resy[i,j]
                - (
                    etas0[i,j] * (
                        (vy1[i,j+1] - vy1[i,j]) / xstp_vy[j]
                        + (vx1[i+1,j] - vx1[i,j]) / ystp_vx[i]
                    )
                    - etas0[i,j-1] * (
                        (vy1[i,j] - vy1[i,j-1]) / xstp_vy[j-1]
                        + (vx1[i+1,j-1] - vx1[i,j-1]) / ystp_vx[i]
                    )
                ) / xstp_b[j-1]
            )
        end
    end
    return nothing
end

function is_valid_node_for_y_stokes(
    i::Int,
    j::Int,
    ynum::Int,
    bintern_zone::Vector{Int}
)::Bool
    return i < ynum + 1 && (j != bintern_zone[5] || i < bintern_zone[6] || i > bintern_zone[7])
end

function is_boundary_node(i::Int, j::Int, xnum::Int, ynum::Int)::Bool
    return i == 1 || i == ynum || j == 1 || j == xnum + 1
end

end # module 