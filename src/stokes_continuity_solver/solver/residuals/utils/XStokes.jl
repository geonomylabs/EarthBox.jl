"""
    XStokes

Compute residual for x-Stokes equation at a given grid location (i,j).

        +----------------------+----------------------+
        |                      |                      |
        |                      |                      |
        |                   vx(i-1,j)                 ------------
        |                      |                      |
        |                      |                      |
        +-----vy(i-1,j)---etas(i-1,j)---vy(i-1,j+1)----    ystpc(i-1)----
        |                      |                      |
        |                      |                      |
    vx(i,j-1)  pr(i-1,j-1)  vx(i,j)     P(i-1,j)   vx(i,j+1)-----   ystp(i-1)
        |    etan(i-1,j-1)     |       etan(i-1,j)    |
        |                      |                      |
        +-------vy(i,j)----etas(i,j)-----vy(i,j+1)-----    ystpc(i)------
        |                      |                      |
        |                      |                      |
        |                   vx(i+1,j)                 -----------
        |                      |                      |
        |                      |                      |
        +----------------------+----------------------+
        |         xstp(j-1)    |      xstp(j)         |
                    |      xstpc(j)        |

Figure 1: X-Stokes equation stencil used to compute x-Stokes residual on
staggered vx grid.
"""
module XStokes

using EarthBox.MathTools
using ..ResidualStructs: StokesResidualsInput

"""
    calculate_residual!(
        residuals_input::StokesResidualsInput,
        i::Int,
        j::Int,
        resx::Matrix{Float64}
    )::Nothing

Compute residual for x-Stokes equation at a given grid location (i,j).

The x-Stokes equation is defined as follows:

    dsxx/dx + dsxy/dy - dP/dx = RX

The residual after solving the system of equations is as follows:

    RX - dsxx/dx - dsxy/dy + dP/dx = residual

See Figure 1 above for details on the stencil used to calculate
derivatives.

# Arguments
- `residuals_input::StokesResidualsInput`: Input parameters for residual calculation
- `i::Int`: Grid index in y direction
- `j::Int`: Grid index in x direction
- `resx::Matrix{Float64}`: Residual array for x-Stokes equation

# Returns
- `Nothing`: Updates resx array in place
"""
function calculate_residual!(
    residuals_input::StokesResidualsInput,
    i::Int,
    j::Int,
    resx::Matrix{Float64}
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
    RX1 = residuals_input.RX1
    vx1 = residuals_input.vx1
    vy1 = residuals_input.vy1
    pr1 = residuals_input.pr1

    if is_valid_node_for_x_stokes(i, j, xnum, ynum, bintern_zone)
        # vx-boundary conditions
        if is_boundary_node(i, j, xnum, ynum)
            resx[i,j] = 0.0
        else
            # dsxx/dx - dP/dx
            resx[i,j] = (
                RX1[i,j]
                - (2.0 * (
                    etan0[i-1,j] * (vx1[i,j+1] - vx1[i,j]) / xstp_b[j]
                    - etan0[i-1,j-1] * (vx1[i,j] - vx1[i,j-1]) / xstp_b[j-1]
                )
                - (pr1[i-1,j] - pr1[i-1,j-1])
                ) / xstp_vy[j]
            )
            # dsxy/dy
            resx[i,j] = (
                resx[i,j]
                - (etas0[i,j] * (
                    (vx1[i+1,j] - vx1[i,j]) / ystp_vx[i]
                    + (vy1[i,j+1] - vy1[i,j]) / xstp_vy[j]
                )
                - etas0[i-1,j] * (
                    (vx1[i,j] - vx1[i-1,j]) / ystp_vx[i-1]
                    + (vy1[i-1,j+1] - vy1[i-1,j]) / xstp_vy[j]
                )
                ) / ystp_b[i-1]
            )
        end
    end
    return nothing
end

function is_valid_node_for_x_stokes(
    i::Int,
    j::Int,
    xnum::Int,
    ynum::Int,
    bintern_zone::Vector{Int}
)::Bool
    return j < xnum + 1 && (j != bintern_zone[1] || i < bintern_zone[2] || i > bintern_zone[3])
end

function is_boundary_node(i::Int, j::Int, xnum::Int, ynum::Int)::Bool
    return i == 1 || i == ynum + 1 || j == 1 || j == xnum
end

end # module 