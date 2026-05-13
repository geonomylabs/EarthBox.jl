module BottomVelocityZNode

import EarthBox.Domain: basic_cell_along_lower_boundary
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
    ynum = build_data.grid.ynum
    ystp = build_data.grid.ystp
    ystpc = build_data.grid.ystpc
    etas_yz = build_data.etas_yz
    bbottomz = build_data.bc.bbottomz
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    if !basic_cell_along_lower_boundary(i, ynum)
        ivzD = ivz + 4
        Lv[inz] = etas_yz[i+1, j, k+1] / ystpc[i+1] / ystp[i]
        Li[inz] = ivz
        Lj[inz] = ivzD
        Lii[inz] = i
        Ljj[inz] = j
        inz = inz + 1
    else
        # vzD is a lower ghost boundary node.
        Lv[inz_c] = (
            Lv[inz_c] +
            bbottomz[j+1, k+1, 2] * etas_yz[i+1, j, k+1] / ystpc[i+1] / ystp[i]
        )
        R[ivz] = (
            R[ivz] -
            bbottomz[j+1, k+1, 1] * etas_yz[i+1, j, k+1] / ystpc[i+1] / ystp[i]
        )
    end
    return inz
end

end # module
