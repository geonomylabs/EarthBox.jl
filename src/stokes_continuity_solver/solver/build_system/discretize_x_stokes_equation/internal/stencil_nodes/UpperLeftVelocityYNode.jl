module UpperLeftVelocityYNode

import EarthBox.Domain: basic_cell_along_upper_boundary
import ....BuildStructs: StokesBuildData
import ....BuildStructs: CellIndices

function coefficients_and_rhs_term(
    inz::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Tuple{Int, Int, Int}
    ivx = cell_indices.ivx
    i = cell_indices.i
    j = cell_indices.j
    ystp = build_data.grid.ystp
    xstpc = build_data.grid.xstpc
    etas = build_data.etas
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
        ivyUL = ivx - 2
        Lv[inz] = etas[i, j+1] / xstpc[j+1] / ystp[i]
        Li[inz] = ivx
        Lj[inz] = ivyUL
        Lii[inz] = i
        Ljj[inz] = j
        # Master non-zero array index for upper left node. This is used when
        # updating coefficients is required.
        inz_UL = inz
        inz = inz + 1
    else
        # Apply boundary conditions since vyUL is a boundary node.
        # The boundary condition equation includes vyLL and,
        # therefore, a term is introduced that requires updating
        # the coefficient for vyLL. A constant is also introduced
        # that must be subtracted from the rhs term of the central node.
        ivyLL = ivx + 1
        Lv[inz] = btopy[j+1, 2] * etas[i, j+1] / xstpc[j+1] / ystp[i]
        Li[inz] = ivx
        Lj[inz] = ivyLL
        Lii[inz] = i
        Ljj[inz] = j
        # Master non-zero array index for lower left node. This is
        # used when updating coefficients is required.
        inz_LL = inz
        inz = inz + 1
        R[ivx] = (
            R[ivx] -
            btopy[j+1, 1] * etas[i, j+1] / xstpc[j+1] / ystp[i]
        )
    end
    return inz, inz_UL, inz_LL
end

end # module 