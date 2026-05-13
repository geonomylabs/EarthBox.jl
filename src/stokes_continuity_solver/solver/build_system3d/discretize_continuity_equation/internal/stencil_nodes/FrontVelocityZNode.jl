module FrontVelocityZNode

import ...Predicates3D: basic_cell_along_front_boundary
import ...Predicates3D: basic_cell_along_back_boundary
import ....BuildStructs: StokesBuildData
import ....BuildStructs: CellIndices

function coefficients_and_rhs_term(
    inz::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Int
    i = cell_indices.i
    j = cell_indices.j
    k = cell_indices.k
    ipr = cell_indices.ipr
    znum = build_data.grid.znum
    zstp = build_data.grid.zstp
    bbackz = build_data.bc.bbackz
    pscale = build_data.pscale
    dshift_to_vxF = build_data.grid.dshift_to_vxF
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    if !basic_cell_along_front_boundary(k)
        ivzF = ipr - 1 - dshift_to_vxF
        Lv[inz] = -pscale / zstp[k]
        Li[inz] = ipr
        Lj[inz] = ivzF
        Lii[inz] = i
        Ljj[inz] = j
        inz_zF = inz
        inz = inz + 1
        # If vzB is the back-boundary ghost, substitute its boundary
        # equation in terms of vzF using bbackz.
        if basic_cell_along_back_boundary(k, znum)
            Lv[inz_zF] = Lv[inz_zF] + bbackz[i+1, j+1, 2] * pscale / zstp[k]
            R[ipr] = R[ipr] - bbackz[i+1, j+1, 1] * pscale / zstp[k]
        end
    end
    return inz
end

end # module
