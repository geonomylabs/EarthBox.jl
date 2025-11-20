"""
    Continuity

Compute residual for continuity equation at a given grid location (i,j).

    +------vy(i,j+1)-------+--------
    |                      |
    |                      |
vx(i+1,j)   pr(i,j)   vx(i+1,j+1) ystp(i)
    |                      |
    |                      |
    +-----vy(i+1,j+1)------+--------
    |        xstp(j)       |

Figure 1: Continuity equation stencil used to calculate continuity 
residual.
"""
module Continuity

using EarthBox.MathTools
using ..ResidualStructs: StokesResidualsInput

"""
    calculate_residual!(
        residuals_input::StokesResidualsInput,
        i::Int,
        j::Int,
        resc::Matrix{Float64}
    )::Nothing

Compute residual for continuity equation at a given grid location.

The continuity equation is defined as follows:
    
    dvx/dx + dvy/dy = RC

The residual after solving the system of equations is as follows:
    
    RC - dvx/dx - dvy/dy = residual

See Figure 1 above for details on the stencil used to calculate 
derivatives.

# Arguments
- `residuals_input::StokesResidualInput`: Input parameters for residual calculation
- `i::Int`: Grid index in y direction
- `j::Int`: Grid index in x direction
- `resc::Array{Float64,2}`: Residual array for continuity equation

# Returns
- `Nothing`: Updates resc array in place
"""
function calculate_residual!(
    residuals_input::StokesResidualsInput,
    i::Int,
    j::Int,
    resc::Matrix{Float64}
)::Nothing
    xnum = residuals_input.xnum
    ynum = residuals_input.ynum
    xstp_b = residuals_input.xstp_b
    ystp_b = residuals_input.ystp_b
    RC1 = residuals_input.RC1
    vx1 = residuals_input.vx1
    vy1 = residuals_input.vy1

    if basic_node_is_not_along_right_or_bottom_edge_of_model(i, j, ynum, xnum)
        resc[i,j] = (
            RC1[i,j] 
            - ((vx1[i+1,j+1] - vx1[i+1,j])/xstp_b[j] 
               + (vy1[i+1,j+1] - vy1[i,j+1])/ystp_b[i])
        )
    end
    return nothing
end

function basic_node_is_not_along_right_or_bottom_edge_of_model(
    i::Int,
    j::Int,
    ynum::Int,
    xnum::Int
)::Bool
    return i < ynum && j < xnum
end

end # module 