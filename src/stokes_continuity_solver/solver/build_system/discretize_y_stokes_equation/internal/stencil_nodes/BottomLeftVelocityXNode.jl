module BottomLeftVelocityXNode

import EarthBox.Domain: basic_cell_along_left_boundary
import ....BuildStructs: StokesBuildData
import ....BuildStructs: CellIndices

function coefficients_and_rhs_term(
    inz::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData,
    drhodx_gy_dt::Float64
)::Tuple{Int,Int,Int}
    i = cell_indices.i
    j = cell_indices.j
    ivy = cell_indices.ivy
    hshift_to_vxR = build_data.grid.hshift_to_vxR
    xstp = build_data.grid.xstp
    ystpc = build_data.grid.ystpc
    etas = build_data.etas
    bleftx = build_data.bc.bleftx
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    if !basic_cell_along_left_boundary(j)
        ivxLL = ivy + 2 - hshift_to_vxR
        Lv[inz] = -etas[i+1, j]/ystpc[i+1]/xstp[j] - drhodx_gy_dt/4.0
        Li[inz] = ivy
        Lj[inz] = ivxLL
        Lii[inz] = i
        Ljj[inz] = j
        inz_LL = inz
        inz = inz + 1
        inz_LR = 0  # Not used in this branch
    else
        # Apply boundary condition equation since vxLL is a ghost node. The
        # boundary condition equation includes vxLR and is described as
        # follows:
        #
        # vxLL = C + D*vxLR
        #
        # where C = bleftx[i+2,1] and D = bleftx[i+2,2]. The coefficient B for
        # vxLL is defined as:
        #
        # B = -etas[i+1,j]/ystpc[i+1]/xstp[j] - drhodx_gy_dt/4.0.
        #
        # When applied to the boundary condition equation, the following is
        # obtained:
        #
        # B*vxLL = B*C + B*D*vxLR
        #
        # Therefore, a term B*D is introduced that requires updating the
        # coefficient for vxLR and constant B*C is introduced that must be
        # subtracted from the rhs term of the central node.
        ivxLR = ivy + 2
        Lv[inz] = (
            bleftx[i+2, 2]*(-etas[i+1, j]/ystpc[i+1]/xstp[j] - drhodx_gy_dt/4.0)
        )
        Li[inz] = ivy
        Lj[inz] = ivxLR
        Lii[inz] = i
        Ljj[inz] = j
        inz_LR = inz
        inz = inz + 1
        inz_LL = 0  # Not used in this branch
        R[ivy] = (
            R[ivy] -
            bleftx[i+2, 1]*(-etas[i+1, j]/ystpc[i+1]/xstp[j] - drhodx_gy_dt/4.0)
        )
    end
    return inz, inz_LL, inz_LR
end

end # module 