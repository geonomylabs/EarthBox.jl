module UpperLeftVelocityXNode

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

    inz_UL = -1
    inz_UR = -1
    if !basic_cell_along_left_boundary(j)
        # ivy - 1 = vx of same cell (i, j, k); back one j-cell is vxUL.
        ivxUL = ivy - 1 - hshift_to_vxR
        Lv[inz] = etas_xy[i+1, j, k] / ystpc[i+1] / xstp[j] - drhodx_gy_dt / 4.0
        Li[inz] = ivy
        Lj[inz] = ivxUL
        Lii[inz] = i
        Ljj[inz] = j
        inz_UL = inz
        inz = inz + 1
    else
        # vxUL is a left ghost boundary node. Substitute the boundary
        # equation that relates vxUL to vxUR:
        #   vxUL = bleftx[i+1, k+1, 1] + bleftx[i+1, k+1, 2]*vxUR
        # producing a coefficient on vxUR and a constant subtracted from
        # the rhs term of the central node.
        ivxUR = ivy - 1
        Lv[inz] = bleftx[i+1, k+1, 2] * (
            etas_xy[i+1, j, k] / ystpc[i+1] / xstp[j] - drhodx_gy_dt / 4.0
        )
        Li[inz] = ivy
        Lj[inz] = ivxUR
        Lii[inz] = i
        Ljj[inz] = j
        inz_UR = inz
        inz = inz + 1
        inz_UL = 0
        R[ivy] = R[ivy] - bleftx[i+1, k+1, 1] * (
            etas_xy[i+1, j, k] / ystpc[i+1] / xstp[j] - drhodx_gy_dt / 4.0
        )
    end
    return inz, inz_UL, inz_UR
end

end # module
