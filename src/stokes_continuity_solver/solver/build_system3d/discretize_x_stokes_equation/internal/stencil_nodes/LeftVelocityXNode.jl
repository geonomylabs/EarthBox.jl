module LeftVelocityXNode

import EarthBox.Domain: basic_cell_along_left_boundary
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
    xstp = build_data.grid.xstp
    xstpc = build_data.grid.xstpc
    etan = build_data.etan
    bleftx = build_data.bc.bleftx
    hshift_to_vxR = build_data.grid.hshift_to_vxR
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    if !basic_cell_along_left_boundary(j)
        # hshift_to_vxR = (ynum-1)*4 in 3D — stride to j+1 neighbour.
        ivxL = ivx - hshift_to_vxR
        Lv[inz] = 2.0 * etan[i, j, k] / xstp[j] / xstpc[j+1]
        Li[inz] = ivx
        Lj[inz] = ivxL
        Lii[inz] = i
        Ljj[inz] = j
        inz = inz + 1
    else
        # Basic node and vxL are along left boundary. vxL is a ghost
        # boundary node so the boundary condition equation
        #     vxL = bleftx[i+1, k+1, 1] + bleftx[i+1, k+1, 2] * vxC
        # is substituted, producing an update to the vxC coefficient and
        # to the rhs term.
        Lv[inz_c] = (
            Lv[inz_c] +
            bleftx[i+1, k+1, 2] * 2.0 * etan[i, j, k] / xstp[j] / xstpc[j+1]
        )
        R[ivx] = (
            R[ivx] -
            bleftx[i+1, k+1, 1] * 2.0 * etan[i, j, k] / xstp[j] / xstpc[j+1]
        )
    end
    return inz
end

end # module
