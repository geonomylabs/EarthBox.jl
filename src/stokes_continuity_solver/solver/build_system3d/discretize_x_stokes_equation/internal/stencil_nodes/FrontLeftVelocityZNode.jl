module FrontLeftVelocityZNode

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
    dshift_to_vxF = build_data.grid.dshift_to_vxF
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    inz_FL = -1
    inz_BL = -1

    if !basic_cell_along_front_boundary(k)
        # vz of cell (i, j, k-1). Within-cell vz offset is +2 from vx;
        # cross-k stride is dshift_to_vxF.
        ivzFL = ivx + 2 - dshift_to_vxF
        Lv[inz] = etas_xz[i, j+1, k] / xstpc[j+1] / zstp[k]
        Li[inz] = ivx
        Lj[inz] = ivzFL
        Lii[inz] = i
        Ljj[inz] = j
        inz_FL = inz
        inz = inz + 1
    else
        # vzFL is at the front boundary. Substitute the boundary equation
        #   vzFL = bfrontz[i+1, j+1, 1] + bfrontz[i+1, j+1, 2]*vzBL
        # producing a coefficient on vzBL (which lives in cell (i, j, k)
        # — the same cell as vxC) and a constant subtracted from the rhs
        # term of the central node.
        ivzBL = ivx + 2
        Lv[inz] = bfrontz[i+1, j+1, 2] * etas_xz[i, j+1, k] / xstpc[j+1] / zstp[k]
        Li[inz] = ivx
        Lj[inz] = ivzBL
        Lii[inz] = i
        Ljj[inz] = j
        inz_BL = inz
        inz = inz + 1
        R[ivx] = (
            R[ivx] -
            bfrontz[i+1, j+1, 1] * etas_xz[i, j+1, k] / xstpc[j+1] / zstp[k]
        )
    end
    return inz, inz_FL, inz_BL
end

end # module
