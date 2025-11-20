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
        # Define the column index for the matrix element utilizing
        # horizontal shift index hshift_to_vxR = (ynum-1)*3
        ivxL = ivx - hshift_to_vxR
        # Define value of matrix element
        Lv[inz] = 2.0 * etan[i, j] / xstp[j] / xstpc[j+1]
        # Define row index for the matrix element
        Li[inz] = ivx
        # Define column index for the matrix element
        Lj[inz] = ivxL
        Lii[inz] = i
        Ljj[inz] = j
        inz = inz + 1
    else
        # Basic node and vxL are along left boundary. For this
        # case apply boundary conditions since vxL is a ghost
        # boundary node.
        #
        # bleftx terms are used in an equation relating the ghost
        # node velocity vxL to internal node velocity vxc.
        #
        #   vxL = bleftx[i+1,1] + bleftx[i+1,2]*vxC
        #
        # where bleftx[i+1,1] is a scalar value representing the
        # prescribed velocity magnitude and bleft[i+1,2] is an
        # integer equal to -1, 0 or 1. For example, if a free slip
        # or no-slip boundary condition is applied bleftx[i+1,1]
        # and bleftx[i+1,2] are equal to 0. However, if a non-zero
        # velocity is applied, bleftx[i+1,1] is non-zero and
        # bleft[i+1,2] is equal to 0.
        #
        # When inserting the above expression into the coefficient
        # for vxL a new term is produced for the coefficient of
        # vxC. Also a new constant term is produced that is
        # subtracted from the right-hand-side term.

        # Adding new term to vxC coefficient
        Lv[inz_c] = (
            Lv[inz_c] +
            bleftx[i+1, 2] * 2.0 * etan[i, j] / xstp[j] / xstpc[j+1]
        )

        # Subtracting new term from rhs term
        R[ivx] = (
            R[ivx] -
            bleftx[i+1, 1] * 2.0 * etan[i, j] / xstp[j] / xstpc[j+1]
        )
    end
    return inz
end

end # module 