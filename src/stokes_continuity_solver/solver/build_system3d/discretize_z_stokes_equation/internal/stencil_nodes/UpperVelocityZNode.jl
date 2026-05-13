module UpperVelocityZNode

import EarthBox.Domain: basic_cell_along_upper_boundary
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
    ystp = build_data.grid.ystp
    ystpc = build_data.grid.ystpc
    etas_yz = build_data.etas_yz
    btopz = build_data.bc.btopz
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    if !basic_cell_along_upper_boundary(i)
        ivzU = ivz - 4
        Lv[inz] = etas_yz[i, j, k+1] / ystpc[i] / ystp[i]
        Li[inz] = ivz
        Lj[inz] = ivzU
        Lii[inz] = i
        Ljj[inz] = j
        inz = inz + 1
    else
        # vzU is an upper ghost boundary node.
        Lv[inz_c] = (
            Lv[inz_c] +
            btopz[j+1, k+1, 2] * etas_yz[i, j, k+1] / ystpc[i] / ystp[i]
        )
        R[ivz] = (
            R[ivz] -
            btopz[j+1, k+1, 1] * etas_yz[i, j, k+1] / ystpc[i] / ystp[i]
        )
    end
    return inz
end

end # module
