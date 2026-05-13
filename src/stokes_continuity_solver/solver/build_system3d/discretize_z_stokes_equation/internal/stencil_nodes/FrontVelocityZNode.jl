module FrontVelocityZNode

import ...Predicates3D: basic_cell_along_front_boundary
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
    zstp = build_data.grid.zstp
    zstpc = build_data.grid.zstpc
    etan = build_data.etan
    bfrontz = build_data.bc.bfrontz
    dshift_to_vxF = build_data.grid.dshift_to_vxF
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    if !basic_cell_along_front_boundary(k)
        # vz of cell (i, j, k-1).
        ivzF = ivz - dshift_to_vxF
        Lv[inz] = 2.0 * etan[i, j, k] / zstp[k] / zstpc[k+1]
        Li[inz] = ivz
        Lj[inz] = ivzF
        Lii[inz] = i
        Ljj[inz] = j
        inz = inz + 1
    else
        # vzF is at the front boundary. Substitute via bfrontz; vzF maps
        # to vzC.
        Lv[inz_c] = (
            Lv[inz_c] +
            bfrontz[i+1, j+1, 2] * 2.0 * etan[i, j, k] / zstp[k] / zstpc[k+1]
        )
        R[ivz] = (
            R[ivz] -
            bfrontz[i+1, j+1, 1] * 2.0 * etan[i, j, k] / zstp[k] / zstpc[k+1]
        )
    end
    return inz
end

end # module
