module FrontLeftVelocityXNode

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

    inz_FL = -1
    inz_FR = -1
    if !basic_cell_along_left_boundary(j)
        # vx of cell (i, j-1, k). Within-cell vx offset from vz is -2;
        # cross-j stride is -hshift_to_vxR.
        ivxFL = ivz - 2 - hshift_to_vxR
        Lv[inz] = etas_xz[i, j, k+1] / xstp[j] / zstpc[k+1]
        Li[inz] = ivz
        Lj[inz] = ivxFL
        Lii[inz] = i
        Ljj[inz] = j
        inz_FL = inz
        inz = inz + 1
    else
        # vxFL is a left ghost boundary node. Substitute via bleftx:
        #   vxFL = bleftx[i+1, k+1, 1] + bleftx[i+1, k+1, 2]*vxFR
        # vxFR is the vx of cell (i, j, k) (same cell as vzC).
        ivxFR = ivz - 2
        Lv[inz] = bleftx[i+1, k+1, 2] * etas_xz[i, j, k+1] / xstp[j] / zstpc[k+1]
        Li[inz] = ivz
        Lj[inz] = ivxFR
        Lii[inz] = i
        Ljj[inz] = j
        inz_FR = inz
        inz = inz + 1
        R[ivz] = (
            R[ivz] -
            bleftx[i+1, k+1, 1] * etas_xz[i, j, k+1] / xstp[j] / zstpc[k+1]
        )
    end
    return inz, inz_FL, inz_FR
end

end # module
