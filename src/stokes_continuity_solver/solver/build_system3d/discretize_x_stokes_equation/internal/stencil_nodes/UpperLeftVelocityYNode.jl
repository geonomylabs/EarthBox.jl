module UpperLeftVelocityYNode

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
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    inz_UL = -1
    inz_LL = -1

    if !basic_cell_along_upper_boundary(i)
        # ivx - 3 = vy of cell (i-1, j, k) in 3D (within-cell stride 4).
        ivyUL = ivx - 3
        Lv[inz] = etas_xy[i, j+1, k] / xstpc[j+1] / ystp[i]
        Li[inz] = ivx
        Lj[inz] = ivyUL
        Lii[inz] = i
        Ljj[inz] = j
        inz_UL = inz
        inz = inz + 1
    else
        # vyUL is a top ghost boundary node. The boundary equation
        #   vyUL = btopy[j+1, k+1, 1] + btopy[j+1, k+1, 2]*vyLL
        # is substituted, producing a coefficient on vyLL and a constant
        # subtracted from the rhs term of the central node.
        ivyLL = ivx + 1
        Lv[inz] = btopy[j+1, k+1, 2] * etas_xy[i, j+1, k] / xstpc[j+1] / ystp[i]
        Li[inz] = ivx
        Lj[inz] = ivyLL
        Lii[inz] = i
        Ljj[inz] = j
        inz_LL = inz
        inz = inz + 1
        R[ivx] = (
            R[ivx] -
            btopy[j+1, k+1, 1] * etas_xy[i, j+1, k] / xstpc[j+1] / ystp[i]
        )
    end
    return inz, inz_UL, inz_LL
end

end # module
