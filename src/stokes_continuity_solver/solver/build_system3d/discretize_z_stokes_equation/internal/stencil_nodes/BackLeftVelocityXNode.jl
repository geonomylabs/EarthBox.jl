module BackLeftVelocityXNode

import EarthBox.Domain: basic_cell_along_left_boundary
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
    hshift_to_vxR = build_data.grid.hshift_to_vxR
    dshift_to_vxF = build_data.grid.dshift_to_vxF
    xstp = build_data.grid.xstp
    zstpc = build_data.grid.zstpc
    etas_xz = build_data.etas_xz
    bleftx = build_data.bc.bleftx
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    inz_BL = -1
    inz_BR = -1
    if !basic_cell_along_left_boundary(j)
        # vx of cell (i, j-1, k+1).
        ivxBL = ivz - 2 - hshift_to_vxR + dshift_to_vxF
        Lv[inz] = -etas_xz[i, j, k+1] / xstp[j] / zstpc[k+1]
        Li[inz] = ivz
        Lj[inz] = ivxBL
        Lii[inz] = i
        Ljj[inz] = j
        inz_BL = inz
        inz = inz + 1
    else
        # vxBL is a left ghost. Substitute via bleftx; vxBR lives in cell
        # (i, j, k+1).
        ivxBR = ivz - 2 + dshift_to_vxF
        Lv[inz] = -bleftx[i+1, k+2, 2] * etas_xz[i, j, k+1] / xstp[j] / zstpc[k+1]
        Li[inz] = ivz
        Lj[inz] = ivxBR
        Lii[inz] = i
        Ljj[inz] = j
        inz_BR = inz
        inz = inz + 1
        R[ivz] = (
            R[ivz] +
            bleftx[i+1, k+2, 1] * etas_xz[i, j, k+1] / xstp[j] / zstpc[k+1]
        )
    end
    return inz, inz_BL, inz_BR
end

end # module
