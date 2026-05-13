module BottomVelocityXNode

import EarthBox.Domain: basic_cell_along_lower_boundary
import ....BuildStructs: StokesBuildData
import ....BuildStructs: CellIndices

function coefficients_and_rhs_term(
    inz::Int,
    inz_c::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Int
    ivx = cell_indices.ivx
    i = cell_indices.i
    j = cell_indices.j
    k = cell_indices.k
    ynum = build_data.grid.ynum
    ystp = build_data.grid.ystp
    ystpc = build_data.grid.ystpc
    etas_xy = build_data.etas_xy
    bbottomx = build_data.bc.bbottomx
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    if !basic_cell_along_lower_boundary(i, ynum)
        ivxD = ivx + 4
        Lv[inz] = etas_xy[i+1, j+1, k] / ystpc[i+1] / ystp[i]
        Li[inz] = ivx
        Lj[inz] = ivxD
        Lii[inz] = i
        Ljj[inz] = j
        inz = inz + 1
    else
        # vxD is a lower ghost boundary node.
        Lv[inz_c] = (
            Lv[inz_c] +
            bbottomx[j+1, k+1, 2] * etas_xy[i+1, j+1, k] / ystpc[i+1] / ystp[i]
        )
        R[ivx] = (
            R[ivx] -
            bbottomx[j+1, k+1, 1] * etas_xy[i+1, j+1, k] / ystpc[i+1] / ystp[i]
        )
    end
    return inz
end

end # module
