module BackVelocityZNode

import ...Predicates3D: basic_cell_in_last_two_back_slabs
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
    ivz = cell_indices.ivz
    znum = build_data.grid.znum
    zstp = build_data.grid.zstp
    zstpc = build_data.grid.zstpc
    etan = build_data.etan
    bbackz = build_data.bc.bbackz
    dshift_to_vxF = build_data.grid.dshift_to_vxF
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    if !basic_cell_in_last_two_back_slabs(k, znum)
        ivzB = ivz + dshift_to_vxF
        Lv[inz] = 2.0 * etan[i, j, k+1] / zstp[k+1] / zstpc[k+1]
        Li[inz] = ivz
        Lj[inz] = ivzB
        Lii[inz] = i
        Ljj[inz] = j
        inz = inz + 1
    else
        # vzB is a back ghost / prescribed-boundary node.
        Lv[inz_c] = (
            Lv[inz_c] +
            bbackz[i+1, j+1, 2] * 2.0 * etan[i, j, k+1] / zstp[k+1] / zstpc[k+1]
        )
        R[ivz] = (
            R[ivz] -
            bbackz[i+1, j+1, 1] * 2.0 * etan[i, j, k+1] / zstp[k+1] / zstpc[k+1]
        )
    end
    return inz
end

end # module
