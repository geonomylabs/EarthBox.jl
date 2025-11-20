module UpperRightVelocityXNode

import EarthBox.Domain: basic_cell_along_right_boundary
import EarthBox.Domain: basic_node_on_left_boundary
import ....BuildStructs: StokesBuildData
import ....BuildStructs: CellIndices

function coefficients_and_rhs_term(
    inz::Int,
    inz_UL::Int,
    inz_UR::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData,
    drhodx_gy_dt::Float64
)::Int
    i = cell_indices.i
    j = cell_indices.j
    ivy = cell_indices.ivy
    xstp = build_data.grid.xstp
    ystpc = build_data.grid.ystpc
    xnum = build_data.grid.xnum
    etas = build_data.etas
    brightx = build_data.bc.brightx
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    if !basic_cell_along_right_boundary(j, xnum)
        ivxUR = ivy - 1
        if !basic_node_on_left_boundary(j)
            Lv[inz] = -etas[i+1, j+1]/ystpc[i+1]/xstp[j] - drhodx_gy_dt/4.0
            Li[inz] = ivy
            Lj[inz] = ivxUR
            Lii[inz] = i
            Ljj[inz] = j
            inz = inz + 1
        else
            # Basic node is at right boundary. Since coefficient
            # was already defined for vxUR we need to update the
            # previous value using the master array index.
            Lv[inz_UR] = (
                Lv[inz_UR] -
                etas[i+1, j+1]/ystpc[i+1]/xstp[j] - drhodx_gy_dt/4.0
            )
        end
    else
        # Apply boundary condition since vxUR is a ghost node. The boundary
        # condition equation includes vxUL and is described as follows:
        #
        # vxUR = C + D*vxUL
        #
        # where C = brightx[i+1,1] and D = brightx[i+1,2]. The coefficient B
        # for vxUR is defined as:
        #
        # B = -etas[i+1,j+1]/ystpc[i+1]/xstp[j] - drhodx_gy_dt/4.0
        #
        # When applied to the boundary condition equation, the following is
        # obtained:
        #
        # B*vxUR = B*C + B*D*vxUL
        #
        # Therefore, a term B*D is introduced that requires updating the
        # coefficient for vxUL and a constant B*C is introduced that must be
        # subtracted from the rhs term of the central node.
        Lv[inz_UL] = (
            Lv[inz_UL] +
            brightx[i+1, 2]*(-etas[i+1, j+1]/ystpc[i+1]/xstp[j] - drhodx_gy_dt/4.0)
        )
        R[ivy] = (
            R[ivy] -
            brightx[i+1, 1]*(-etas[i+1, j+1]/ystpc[i+1]/xstp[j] - drhodx_gy_dt/4.0)
        )
    end
    return inz
end

end # module 