module BackLeftVelocityZNode

import ...Predicates3D: basic_cell_along_front_boundary
import ...Predicates3D: basic_cell_along_back_boundary
import ....BuildStructs: StokesBuildData
import ....BuildStructs: CellIndices

function coefficients_and_rhs_term(
    inz::Int,
    inz_FL::Int,
    inz_BL::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Int
    ivx = cell_indices.ivx
    i = cell_indices.i
    j = cell_indices.j
    k = cell_indices.k
    znum = build_data.grid.znum
    xstpc = build_data.grid.xstpc
    zstp = build_data.grid.zstp
    etas_xz = build_data.etas_xz
    bbackz = build_data.bc.bbackz
    dshift_to_vxF = build_data.grid.dshift_to_vxF
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    if !basic_cell_along_back_boundary(k, znum)
        # vz of cell (i, j, k) — same cell as vxC.
        ivzBL = ivx + 2
        if !basic_cell_along_front_boundary(k)
            Lv[inz] = -etas_xz[i, j+1, k+1] / xstpc[j+1] / zstp[k]
            Li[inz] = ivx
            Lj[inz] = ivzBL
            Lii[inz] = i
            Ljj[inz] = j
            inz = inz + 1
        else
            # Basic node is at front boundary; vzBL coefficient was already
            # defined by FrontLeftVelocityZNode at inz_BL. Update it.
            Lv[inz_BL] = (
                Lv[inz_BL] -
                etas_xz[i, j+1, k+1] / xstpc[j+1] / zstp[k]
            )
        end
    else
        # vzBL is a back ghost boundary node. Substitute via bbackz:
        #   vzBL = C + D*vzFL
        # which updates the vzFL coefficient at inz_FL and the rhs term.
        Lv[inz_FL] = (
            Lv[inz_FL] -
            bbackz[i+1, j+1, 2] * etas_xz[i, j+1, k+1] / xstpc[j+1] / zstp[k]
        )
        R[ivx] = (
            R[ivx] +
            bbackz[i+1, j+1, 1] * etas_xz[i, j+1, k+1] / xstpc[j+1] / zstp[k]
        )
    end
    return inz
end

end # module
