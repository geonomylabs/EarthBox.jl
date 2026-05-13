module BackVelocityYNode

import ...Predicates3D: basic_cell_along_back_boundary
import ....BuildStructs: StokesBuildData
import ....BuildStructs: CellIndices

function coefficients_and_rhs_term(
    inz::Int,
    inz_c::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Int
    i = cell_indices.i
    j = cell_indices.j
    k = cell_indices.k
    ivy = cell_indices.ivy
    znum = build_data.grid.znum
    zstp = build_data.grid.zstp
    zstpc = build_data.grid.zstpc
    etas_yz = build_data.etas_yz
    bbacky = build_data.bc.bbacky
    dshift_to_vxF = build_data.grid.dshift_to_vxF
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    if !basic_cell_along_back_boundary(k, znum)
        ivyB = ivy + dshift_to_vxF
        Lv[inz] = etas_yz[i+1, j, k+1] / zstpc[k+1] / zstp[k]
        Li[inz] = ivy
        Lj[inz] = ivyB
        Lii[inz] = i
        Ljj[inz] = j
        inz = inz + 1
    else
        # vyB is a back ghost boundary node.
        Lv[inz_c] = (
            Lv[inz_c] +
            bbacky[i+1, j+1, 2] * etas_yz[i+1, j, k+1] / zstpc[k+1] / zstp[k]
        )
        R[ivy] = (
            R[ivy] -
            bbacky[i+1, j+1, 1] * etas_yz[i+1, j, k+1] / zstpc[k+1] / zstp[k]
        )
    end
    return inz
end

end # module
