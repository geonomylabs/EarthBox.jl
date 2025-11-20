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
        ivyU = ivy - 3
        Lv[inz] = 2.0 * etan[i, j] / ystp[i] / ystpc[i+1]
        Li[inz] = ivy
        Lj[inz] = ivyU
        Lii[inz] = i
        Ljj[inz] = j
        inz = inz + 1
    else
        # Apply boundary condition since vyU is a ghost boundary
        # node. This leads to new terms for the coefficient for
        # vyC and the rhs term.
        Lv[inz_c] = (
            Lv[inz_c] +
            btopy[j+1, 2] * 2.0 * etan[i, j] / ystp[i] / ystpc[i+1]
        )
        R[ivy] = (
            R[ivy] -
            btopy[j+1, 1] * 2.0 * etan[i, j] / ystp[i] / ystpc[i+1]
        )
    end
    return inz
end

end # module 