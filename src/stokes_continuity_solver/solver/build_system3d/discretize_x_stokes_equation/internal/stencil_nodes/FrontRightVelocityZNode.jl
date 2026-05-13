module FrontRightVelocityZNode

import ...Predicates3D: basic_cell_along_front_boundary
import ....BuildStructs: StokesBuildData
import ....BuildStructs: CellIndices

function coefficients_and_rhs_term(
    inz::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Tuple{Int,Int,Int}
    ivx = cell_indices.ivx
    i = cell_indices.i
    j = cell_indices.j
    k = cell_indices.k
    zstp = build_data.grid.zstp
    xstpc = build_data.grid.xstpc
    etas_xz = build_data.etas_xz
    bfrontz = build_data.bc.bfrontz
    hshift_to_vxR = build_data.grid.hshift_to_vxR
    dshift_to_vxF = build_data.grid.dshift_to_vxF
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    inz_FR = -1
    inz_BR = -1

    if !basic_cell_along_front_boundary(k)
        # vz of cell (i, j+1, k-1).
        ivzFR = ivx + 2 + hshift_to_vxR - dshift_to_vxF
        Lv[inz] = -etas_xz[i, j+1, k] / xstpc[j+1] / zstp[k]
        Li[inz] = ivx
        Lj[inz] = ivzFR
        Lii[inz] = i
        Ljj[inz] = j
        inz_FR = inz
        inz = inz + 1
    else
        # vzFR is at the front boundary. Substitute its boundary equation
        # in terms of vzBR (which lives in cell (i, j+1, k)).
        ivzBR = ivx + 2 + hshift_to_vxR
        Lv[inz] = -bfrontz[i+1, j+2, 2] * etas_xz[i, j+1, k] / xstpc[j+1] / zstp[k]
        Li[inz] = ivx
        Lj[inz] = ivzBR
        Lii[inz] = i
        Ljj[inz] = j
        inz_BR = inz
        inz = inz + 1
        R[ivx] = (
            R[ivx] +
            bfrontz[i+1, j+2, 1] * etas_xz[i, j+1, k] / xstpc[j+1] / zstp[k]
        )
    end
    return inz, inz_FR, inz_BR
end

end # module
