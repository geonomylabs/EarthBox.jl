module BackRightVelocityZNode

import ...Predicates3D: basic_cell_along_front_boundary
import ...Predicates3D: basic_cell_along_back_boundary
import ....BuildStructs: StokesBuildData
import ....BuildStructs: CellIndices

function coefficients_and_rhs_term(
    inz::Int,
    inz_FR::Int,
    inz_BR::Int,
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
    hshift_to_vxR = build_data.grid.hshift_to_vxR
    dshift_to_vxF = build_data.grid.dshift_to_vxF
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    if !basic_cell_along_back_boundary(k, znum)
        # vz of cell (i, j+1, k).
        ivzBR = ivx + 2 + hshift_to_vxR
        if !basic_cell_along_front_boundary(k)
            Lv[inz] = etas_xz[i, j+1, k+1] / xstpc[j+1] / zstp[k]
            Li[inz] = ivx
            Lj[inz] = ivzBR
            Lii[inz] = i
            Ljj[inz] = j
            inz = inz + 1
        else
            # Basic node is at front boundary; vzBR coefficient was already
            # defined by FrontRightVelocityZNode at inz_BR. Update it.
            Lv[inz_BR] = (
                Lv[inz_BR] +
                etas_xz[i, j+1, k+1] / xstpc[j+1] / zstp[k]
            )
        end
    else
        # vzBR is a back ghost boundary node. Substitute via bbackz:
        #   vzBR = C + D*vzFR
        Lv[inz_FR] = (
            Lv[inz_FR] +
            bbackz[i+1, j+2, 2] * etas_xz[i, j+1, k+1] / xstpc[j+1] / zstp[k]
        )
        R[ivx] = (
            R[ivx] -
            bbackz[i+1, j+2, 1] * etas_xz[i, j+1, k+1] / xstpc[j+1] / zstp[k]
        )
    end
    return inz
end

end # module
