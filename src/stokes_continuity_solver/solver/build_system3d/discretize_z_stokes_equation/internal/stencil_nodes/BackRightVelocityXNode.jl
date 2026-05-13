module BackRightVelocityXNode

import EarthBox.Domain: basic_cell_along_left_boundary
import EarthBox.Domain: basic_cell_along_right_boundary
import ....BuildStructs: StokesBuildData
import ....BuildStructs: CellIndices

function coefficients_and_rhs_term(
    inz::Int,
    inz_BL::Int,
    inz_BR::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Int
    i = cell_indices.i
    j = cell_indices.j
    k = cell_indices.k
    ivz = cell_indices.ivz
    xnum = build_data.grid.xnum
    dshift_to_vxF = build_data.grid.dshift_to_vxF
    xstp = build_data.grid.xstp
    zstpc = build_data.grid.zstpc
    etas_xz = build_data.etas_xz
    brightx = build_data.bc.brightx
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    if !basic_cell_along_right_boundary(j, xnum)
        ivxBR = ivz - 2 + dshift_to_vxF
        if !basic_cell_along_left_boundary(j)
            Lv[inz] = etas_xz[i, j+1, k+1] / xstp[j] / zstpc[k+1]
            Li[inz] = ivz
            Lj[inz] = ivxBR
            Lii[inz] = i
            Ljj[inz] = j
            inz = inz + 1
        else
            # Basic cell is at left boundary; vxBR coefficient was already
            # defined by BackLeftVelocityXNode at inz_BR. Update it.
            Lv[inz_BR] = (
                Lv[inz_BR] +
                etas_xz[i, j+1, k+1] / xstp[j] / zstpc[k+1]
            )
        end
    else
        # vxBR is a right ghost. Substitute via brightx.
        Lv[inz_BL] = (
            Lv[inz_BL] +
            brightx[i+1, k+2, 2] * etas_xz[i, j+1, k+1] / xstp[j] / zstpc[k+1]
        )
        R[ivz] = (
            R[ivz] -
            brightx[i+1, k+2, 1] * etas_xz[i, j+1, k+1] / xstp[j] / zstpc[k+1]
        )
    end
    return inz
end

end # module
