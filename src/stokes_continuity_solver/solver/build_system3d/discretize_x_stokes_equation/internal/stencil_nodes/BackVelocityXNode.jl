module BackVelocityXNode

import ...Predicates3D: basic_cell_along_back_boundary
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
    znum = build_data.grid.znum
    zstp = build_data.grid.zstp
    zstpc = build_data.grid.zstpc
    etas_xz = build_data.etas_xz
    bbackx = build_data.bc.bbackx
    dshift_to_vxF = build_data.grid.dshift_to_vxF
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    if !basic_cell_along_back_boundary(k, znum)
        ivxB = ivx + dshift_to_vxF
        Lv[inz] = etas_xz[i, j+1, k+1] / zstpc[k+1] / zstp[k]
        Li[inz] = ivx
        Lj[inz] = ivxB
        Lii[inz] = i
        Ljj[inz] = j
        inz = inz + 1
    else
        # vxB is a back ghost boundary node.
        Lv[inz_c] = (
            Lv[inz_c] +
            bbackx[i+1, j+1, 2] * etas_xz[i, j+1, k+1] / zstpc[k+1] / zstp[k]
        )
        R[ivx] = (
            R[ivx] -
            bbackx[i+1, j+1, 1] * etas_xz[i, j+1, k+1] / zstpc[k+1] / zstp[k]
        )
    end
    return inz
end

end # module
