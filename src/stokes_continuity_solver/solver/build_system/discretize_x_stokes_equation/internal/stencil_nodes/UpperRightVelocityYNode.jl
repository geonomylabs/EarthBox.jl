module UpperRightVelocityYNode

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
        ivyUR = ivx - 2 + hshift_to_vxR
        Lv[inz] = -etas[i, j+1] / xstpc[j+1] / ystp[i]
        Li[inz] = ivx
        Lj[inz] = ivyUR
        Lii[inz] = i
        Ljj[inz] = j
        inz_UR = inz
        inz = inz + 1
    else
        # Apply boundary conditions since vyUR is a boundary ghost
        # node. The boundary condition equation includes vyLR and,
        # therefore, a term is introduced that requires updating
        # the coefficient for vyLR. A constant is also introduced
        # that must be subtracted from the rhs term of the central node.
        ivyLR = ivx + 1 + hshift_to_vxR
        Lv[inz] = -btopy[j+2, 2] * etas[i, j+1] / xstpc[j+1] / ystp[i]
        Li[inz] = ivx
        Lj[inz] = ivyLR
        Lii[inz] = i
        Ljj[inz] = j
        # Master non-zero array index for lower right node. This
        # is used when updating coefficients is required.
        inz_LR = inz
        inz = inz + 1
        R[ivx] = (
            R[ivx] +
            btopy[j+2, 1] * etas[i, j+1] / xstpc[j+1] / ystp[i]
        )
    end
    return inz, inz_UR, inz_LR
end

end # module 