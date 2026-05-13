module UpperVelocityYNode

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
    ivy = cell_indices.ivy
    ystp = build_data.grid.ystp
    ystpc = build_data.grid.ystpc
    etan = build_data.etan
    btopy = build_data.bc.btopy
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    if !basic_cell_along_upper_boundary(i)
        # vy of cell (i-1, j, k); within-cell stride is 4.
        ivyU = ivy - 4
        Lv[inz] = 2.0 * etan[i, j, k] / ystp[i] / ystpc[i+1]
        Li[inz] = ivy
        Lj[inz] = ivyU
        Lii[inz] = i
        Ljj[inz] = j
        inz = inz + 1
    else
        # vyU is an upper ghost boundary node.
        Lv[inz_c] = (
            Lv[inz_c] +
            btopy[j+1, k+1, 2] * 2.0 * etan[i, j, k] / ystp[i] / ystpc[i+1]
        )
        R[ivy] = (
            R[ivy] -
            btopy[j+1, k+1, 1] * 2.0 * etan[i, j, k] / ystp[i] / ystpc[i+1]
        )
    end
    return inz
end

end # module
