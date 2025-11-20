module BottomRightVelocityXNode

import EarthBox.Domain: basic_cell_along_left_boundary
import EarthBox.Domain: basic_cell_along_right_boundary
import ....BuildStructs: StokesBuildData
import ....BuildStructs: CellIndices

function coefficients_and_rhs_term(
    inz::Int,
    inz_LL::Int,
    inz_LR::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData,
    drhodx_gy_dt::Float64
)::Int
    i = cell_indices.i
    j = cell_indices.j
    ivy = cell_indices.ivy
    xnum = build_data.grid.xnum
    xstp = build_data.grid.xstp
    ystpc = build_data.grid.ystpc
    etas = build_data.etas
    brightx = build_data.bc.brightx
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    if !basic_cell_along_right_boundary(j, xnum)
        ivxLR = ivy + 2
        if !basic_cell_along_left_boundary(j)
            Lv[inz] = etas[i+1, j+1]/ystpc[i+1]/xstp[j] - drhodx_gy_dt/4.0
            Li[inz] = ivy
            Lj[inz] = ivxLR
            Lii[inz] = i
            Ljj[inz] = j
            inz = inz + 1
        else
            # Basic node is at right boundary. Since coefficient
            # was already defined for vxLR we need to update the
            # previous value using the master array index.
            Lv[inz_LR] = (
                Lv[inz_LR] +
                etas[i+1, j+1]/ystpc[i+1]/xstp[j] - drhodx_gy_dt/4.0
            )
        end
    else
        # Apply boundary condition since vxLR is a ghost node. The boundary
        # condition equation includes vxLL and is described as follows:
        #
        # vxLR = C + D*vxLL
        #
        # where C = brightx[i+1,1] and D = brightx[i+1,2]. The coefficient B
        # for vxLR is defined as:
        #
        # B = etas[i+1,j+1]/ystpc[i+1]/xstp[j] - drhodx_gy_dt/4.0
        #
        # When applied to the boundary condition equation, the following is
        # obtained:
        #
        # B*vxLR = B*C + B*D*vxLL
        #
        # Therefore, a term is introduced that requires updating the
        # coefficient for vxLL. A constant is also introduced that must be
        # added to the rhs term of the central node.
        Lv[inz_LL] = (
            Lv[inz_LL] +
            brightx[i+2, 2]*(etas[i+1, j+1]/ystpc[i+1]/xstp[j] - drhodx_gy_dt/4.0)
        )
        R[ivy] = (
            R[ivy] -
            brightx[i+2, 1]*(etas[i+1, j+1]/ystpc[i+1]/xstp[j] - drhodx_gy_dt/4.0)
        )
    end
    return inz
end

end # module 