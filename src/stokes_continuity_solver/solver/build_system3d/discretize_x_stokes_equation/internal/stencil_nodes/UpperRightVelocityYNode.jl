module UpperRightVelocityYNode

import EarthBox.Domain: basic_cell_along_upper_boundary
import ....BuildStructs: StokesBuildData
import ....BuildStructs: CellIndices

function coefficients_and_rhs_term(
    inz::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Tuple{Int,Int,Int}
    ivx = cell_indices.ivx
    i = cell_indices.i
    j = cell_indices.j
    k = cell_indices.k
    ystp = build_data.grid.ystp
    xstpc = build_data.grid.xstpc
    etas_xy = build_data.etas_xy
    btopy = build_data.bc.btopy
    hshift_to_vxR = build_data.grid.hshift_to_vxR
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    inz_UR = -1
    inz_LR = -1

    if !basic_cell_along_upper_boundary(i)
        ivyUR = ivx - 3 + hshift_to_vxR
        Lv[inz] = -etas_xy[i, j+1, k] / xstpc[j+1] / ystp[i]
        Li[inz] = ivx
        Lj[inz] = ivyUR
        Lii[inz] = i
        Ljj[inz] = j
        inz_UR = inz
        inz = inz + 1
    else
        # vyUR is a top ghost boundary node. Substitute the boundary
        # equation that relates vyUR to vyLR.
        ivyLR = ivx + 1 + hshift_to_vxR
        Lv[inz] = -btopy[j+2, k+1, 2] * etas_xy[i, j+1, k] / xstpc[j+1] / ystp[i]
        Li[inz] = ivx
        Lj[inz] = ivyLR
        Lii[inz] = i
        Ljj[inz] = j
        inz_LR = inz
        inz = inz + 1
        R[ivx] = (
            R[ivx] +
            btopy[j+2, k+1, 1] * etas_xy[i, j+1, k] / xstpc[j+1] / ystp[i]
        )
    end
    return inz, inz_UR, inz_LR
end

end # module
