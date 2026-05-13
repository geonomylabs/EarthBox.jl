module LeftVelocityZNode

import EarthBox.Domain: basic_cell_along_left_boundary
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
    xstp = build_data.grid.xstp
    xstpc = build_data.grid.xstpc
    etas_xz = build_data.etas_xz
    bleftz = build_data.bc.bleftz
    hshift_to_vxR = build_data.grid.hshift_to_vxR
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    if !basic_cell_along_left_boundary(j)
        ivzL = ivz - hshift_to_vxR
        Lv[inz] = etas_xz[i, j, k+1] / xstpc[j] / xstp[j]
        Li[inz] = ivz
        Lj[inz] = ivzL
        Lii[inz] = i
        Ljj[inz] = j
        inz = inz + 1
    else
        # vzL is a left ghost boundary node.
        Lv[inz_c] = (
            Lv[inz_c] +
            bleftz[i+1, k+1, 2] * etas_xz[i, j, k+1] / xstpc[j] / xstp[j]
        )
        R[ivz] = (
            R[ivz] -
            bleftz[i+1, k+1, 1] * etas_xz[i, j, k+1] / xstpc[j] / xstp[j]
        )
    end
    return inz
end

end # module
