module FrontUpperVelocityYNode

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
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    inz_FU = -1
    inz_FD = -1

    if !basic_cell_along_upper_boundary(i)
        # vy of cell (i-1, j, k). Within-cell vy offset from vz is -1;
        # cross-i stride is -4. Net: ivz - 1 - 4 = ivz - 5.
        ivyFU = ivz - 5
        Lv[inz] = etas_yz[i, j, k+1] / ystp[i] / zstpc[k+1]
        Li[inz] = ivz
        Lj[inz] = ivyFU
        Lii[inz] = i
        Ljj[inz] = j
        inz_FU = inz
        inz = inz + 1
    else
        # vyFU is at the top boundary. Substitute via btopy:
        #   vyFU = btopy[j+1, k+1, 1] + btopy[j+1, k+1, 2]*vyFD
        # vyFD is the vy of cell (i, j, k) (same cell as vzC).
        ivyFD = ivz - 1
        Lv[inz] = btopy[j+1, k+1, 2] * etas_yz[i, j, k+1] / ystp[i] / zstpc[k+1]
        Li[inz] = ivz
        Lj[inz] = ivyFD
        Lii[inz] = i
        Ljj[inz] = j
        inz_FD = inz
        inz = inz + 1
        R[ivz] = (
            R[ivz] -
            btopy[j+1, k+1, 1] * etas_yz[i, j, k+1] / ystp[i] / zstpc[k+1]
        )
    end
    return inz, inz_FU, inz_FD
end

end # module
