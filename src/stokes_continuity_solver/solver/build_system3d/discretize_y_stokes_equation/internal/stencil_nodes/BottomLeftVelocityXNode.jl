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
    k = cell_indices.k
    ivy = cell_indices.ivy
    hshift_to_vxR = build_data.grid.hshift_to_vxR
    xstp = build_data.grid.xstp
    ystpc = build_data.grid.ystpc
    etas_xy = build_data.etas_xy
    bleftx = build_data.bc.bleftx
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    if !basic_cell_along_left_boundary(j)
        # vx of cell (i+1, j-1, k). Within-cell vx offset is -1 from ivy.
        # Same-axis stride between i and i+1 is +4. Cross-axis (j-1) is
        # -hshift_to_vxR. Net: ivy - 1 + 4 - hshift = ivy + 3 - hshift.
        ivxLL = ivy + 3 - hshift_to_vxR
        Lv[inz] = -etas_xy[i+1, j, k] / ystpc[i+1] / xstp[j] - drhodx_gy_dt / 4.0
        Li[inz] = ivy
        Lj[inz] = ivxLL
        Lii[inz] = i
        Ljj[inz] = j
        inz_LL = inz
        inz = inz + 1
        inz_LR = 0
    else
        # vxLL is a left ghost. Substitute via bleftx: vxLL = C + D*vxLR.
        ivxLR = ivy + 3
        Lv[inz] = bleftx[i+2, k+1, 2] * (
            -etas_xy[i+1, j, k] / ystpc[i+1] / xstp[j] - drhodx_gy_dt / 4.0
        )
        Li[inz] = ivy
        Lj[inz] = ivxLR
        Lii[inz] = i
        Ljj[inz] = j
        inz_LR = inz
        inz = inz + 1
        inz_LL = 0
        R[ivy] = R[ivy] - bleftx[i+2, k+1, 1] * (
            -etas_xy[i+1, j, k] / ystpc[i+1] / xstp[j] - drhodx_gy_dt / 4.0
        )
    end
    return inz, inz_LL, inz_LR
end

end # module
