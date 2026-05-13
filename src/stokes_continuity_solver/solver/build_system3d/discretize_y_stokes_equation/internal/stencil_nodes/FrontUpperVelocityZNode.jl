module FrontUpperVelocityZNode

import ...Predicates3D: basic_cell_along_front_boundary
import ....BuildStructs: StokesBuildData
import ....BuildStructs: CellIndices

function coefficients_and_rhs_term(
    inz::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData,
    drhodz_gy_dt::Float64
)::Tuple{Int,Int,Int}
    i = cell_indices.i
    j = cell_indices.j
    k = cell_indices.k
    ivy = cell_indices.ivy
    dshift_to_vxF = build_data.grid.dshift_to_vxF
    ystpc = build_data.grid.ystpc
    zstp = build_data.grid.zstp
    etas_yz = build_data.etas_yz
    bfrontz = build_data.bc.bfrontz
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    inz_FU = -1
    inz_BU = -1
    if !basic_cell_along_front_boundary(k)
        # vz of cell (i, j, k-1). Within-cell vz offset from vy is +1.
        ivzFU = ivy + 1 - dshift_to_vxF
        Lv[inz] = etas_yz[i+1, j, k] / ystpc[i+1] / zstp[k] - drhodz_gy_dt / 4.0
        Li[inz] = ivy
        Lj[inz] = ivzFU
        Lii[inz] = i
        Ljj[inz] = j
        inz_FU = inz
        inz = inz + 1
    else
        # vzFU is at the front boundary. Substitute via bfrontz:
        #   vzFU = bfrontz[i+1, j+1, 1] + bfrontz[i+1, j+1, 2]*vzBU
        # vzBU lives in cell (i, j, k) — the same cell as vyC.
        ivzBU = ivy + 1
        Lv[inz] = bfrontz[i+1, j+1, 2] * (
            etas_yz[i+1, j, k] / ystpc[i+1] / zstp[k] - drhodz_gy_dt / 4.0
        )
        Li[inz] = ivy
        Lj[inz] = ivzBU
        Lii[inz] = i
        Ljj[inz] = j
        inz_BU = inz
        inz = inz + 1
        inz_FU = 0
        R[ivy] = R[ivy] - bfrontz[i+1, j+1, 1] * (
            etas_yz[i+1, j, k] / ystpc[i+1] / zstp[k] - drhodz_gy_dt / 4.0
        )
    end
    return inz, inz_FU, inz_BU
end

end # module
