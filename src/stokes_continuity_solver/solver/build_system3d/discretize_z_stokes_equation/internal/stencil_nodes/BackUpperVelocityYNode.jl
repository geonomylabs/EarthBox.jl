module BackUpperVelocityYNode

import EarthBox.Domain: basic_cell_along_upper_boundary
import ....BuildStructs: StokesBuildData
import ....BuildStructs: CellIndices

function coefficients_and_rhs_term(
    inz::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Tuple{Int,Int,Int}
    i = cell_indices.i
    j = cell_indices.j
    k = cell_indices.k
    ivz = cell_indices.ivz
    ystp = build_data.grid.ystp
    zstpc = build_data.grid.zstpc
    etas_yz = build_data.etas_yz
    btopy = build_data.bc.btopy
    dshift_to_vxF = build_data.grid.dshift_to_vxF
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    inz_BU = -1
    inz_BD = -1

    if !basic_cell_along_upper_boundary(i)
        # vy of cell (i-1, j, k+1).
        ivyBU = ivz - 5 + dshift_to_vxF
        Lv[inz] = -etas_yz[i, j, k+1] / ystp[i] / zstpc[k+1]
        Li[inz] = ivz
        Lj[inz] = ivyBU
        Lii[inz] = i
        Ljj[inz] = j
        inz_BU = inz
        inz = inz + 1
    else
        # vyBU at top boundary. Substitute via btopy. vyBD lives in cell
        # (i, j, k+1).
        ivyBD = ivz - 1 + dshift_to_vxF
        Lv[inz] = -btopy[j+1, k+2, 2] * etas_yz[i, j, k+1] / ystp[i] / zstpc[k+1]
        Li[inz] = ivz
        Lj[inz] = ivyBD
        Lii[inz] = i
        Ljj[inz] = j
        inz_BD = inz
        inz = inz + 1
        R[ivz] = (
            R[ivz] +
            btopy[j+1, k+2, 1] * etas_yz[i, j, k+1] / ystp[i] / zstpc[k+1]
        )
    end
    return inz, inz_BU, inz_BD
end

end # module
