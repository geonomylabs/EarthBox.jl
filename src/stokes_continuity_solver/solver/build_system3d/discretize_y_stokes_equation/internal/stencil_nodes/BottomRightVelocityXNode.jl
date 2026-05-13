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
    k = cell_indices.k
    ivy = cell_indices.ivy
    xnum = build_data.grid.xnum
    xstp = build_data.grid.xstp
    ystpc = build_data.grid.ystpc
    etas_xy = build_data.etas_xy
    brightx = build_data.bc.brightx
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    if !basic_cell_along_right_boundary(j, xnum)
        ivxLR = ivy + 3
        if !basic_cell_along_left_boundary(j)
            Lv[inz] = etas_xy[i+1, j+1, k] / ystpc[i+1] / xstp[j] - drhodx_gy_dt / 4.0
            Li[inz] = ivy
            Lj[inz] = ivxLR
            Lii[inz] = i
            Ljj[inz] = j
            inz = inz + 1
        else
            # Basic node is at left boundary; vxLR coefficient was already
            # defined by BottomLeftVelocityXNode at inz_LR. Update it.
            Lv[inz_LR] = (
                Lv[inz_LR] +
                etas_xy[i+1, j+1, k] / ystpc[i+1] / xstp[j] - drhodx_gy_dt / 4.0
            )
        end
    else
        # vxLR is a right ghost boundary node. Substitute via brightx.
        Lv[inz_LL] = Lv[inz_LL] + brightx[i+2, k+1, 2] * (
            etas_xy[i+1, j+1, k] / ystpc[i+1] / xstp[j] - drhodx_gy_dt / 4.0
        )
        R[ivy] = R[ivy] - brightx[i+2, k+1, 1] * (
            etas_xy[i+1, j+1, k] / ystpc[i+1] / xstp[j] - drhodx_gy_dt / 4.0
        )
    end
    return inz
end

end # module
