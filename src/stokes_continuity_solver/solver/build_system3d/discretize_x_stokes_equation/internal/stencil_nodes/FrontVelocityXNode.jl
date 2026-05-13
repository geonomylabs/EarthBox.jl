module FrontVelocityXNode

import ...Predicates3D: basic_cell_along_front_boundary
import ....BuildStructs: StokesBuildData
import ....BuildStructs: CellIndices

function coefficients_and_rhs_term(
    inz::Int,
    inz_c::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Int
    ivx = cell_indices.ivx
    i = cell_indices.i
    j = cell_indices.j
    k = cell_indices.k
    zstp = build_data.grid.zstp
    zstpc = build_data.grid.zstpc
    etas_xz = build_data.etas_xz
    bfrontx = build_data.bc.bfrontx
    dshift_to_vxF = build_data.grid.dshift_to_vxF
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    if !basic_cell_along_front_boundary(k)
        # dshift_to_vxF = (xnum-1)*(ynum-1)*4 — stride to k-1 neighbour.
        ivxF = ivx - dshift_to_vxF
        Lv[inz] = etas_xz[i, j+1, k] / zstpc[k] / zstp[k]
        Li[inz] = ivx
        Lj[inz] = ivxF
        Lii[inz] = i
        Ljj[inz] = j
        inz = inz + 1
    else
        # vxF is a front ghost boundary node.
        Lv[inz_c] = (
            Lv[inz_c] +
            bfrontx[i+1, j+1, 2] * etas_xz[i, j+1, k] / zstpc[k] / zstp[k]
        )
        R[ivx] = (
            R[ivx] -
            bfrontx[i+1, j+1, 1] * etas_xz[i, j+1, k] / zstpc[k] / zstp[k]
        )
    end
    return inz
end

end # module
