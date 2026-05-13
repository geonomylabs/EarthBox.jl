module RightVelocityZNode

import EarthBox.Domain: basic_cell_along_right_boundary
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
    xnum = build_data.grid.xnum
    xstp = build_data.grid.xstp
    xstpc = build_data.grid.xstpc
    etas_xz = build_data.etas_xz
    brightz = build_data.bc.brightz
    hshift_to_vxR = build_data.grid.hshift_to_vxR
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    if !basic_cell_along_right_boundary(j, xnum)
        ivzR = ivz + hshift_to_vxR
        Lv[inz] = etas_xz[i, j+1, k+1] / xstpc[j+1] / xstp[j]
        Li[inz] = ivz
        Lj[inz] = ivzR
        Lii[inz] = i
        Ljj[inz] = j
        inz = inz + 1
    else
        # vzR is a right ghost boundary node.
        Lv[inz_c] = (
            Lv[inz_c] +
            brightz[i+1, k+1, 2] * etas_xz[i, j+1, k+1] / xstpc[j+1] / xstp[j]
        )
        R[ivz] = (
            R[ivz] -
            brightz[i+1, k+1, 1] * etas_xz[i, j+1, k+1] / xstpc[j+1] / xstp[j]
        )
    end
    return inz
end

end # module
