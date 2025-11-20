module LeftVelocityXNode

import EarthBox.Domain: basic_cell_along_left_boundary
import EarthBox.Domain: basic_cell_along_right_boundary
import ....BuildStructs: StokesBuildData
import ....BuildStructs: CellIndices

function coefficients_and_rhs_term(
    inz::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Int
    i = cell_indices.i
    j = cell_indices.j
    ipr = cell_indices.ipr
    xnum = build_data.grid.xnum
    xstp = build_data.grid.xstp
    brightx = build_data.bc.brightx
    pscale = build_data.pscale
    hshift_to_vxR = build_data.grid.hshift_to_vxR
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    if !basic_cell_along_left_boundary(j)
        ivxL = ipr - 2 - hshift_to_vxR
        Lv[inz] = -pscale / xstp[j]
        Li[inz] = ipr
        Lj[inz] = ivxL
        Lii[inz] = i
        Ljj[inz] = j
        inz_xL = inz
        inz = inz + 1
        # Add boundary condition for the right Vx node
        if basic_cell_along_right_boundary(j, xnum)
            Lv[inz_xL] = Lv[inz_xL] + brightx[i+1, 2] * pscale / xstp[j]
            R[ipr] = R[ipr] - brightx[i+1, 1] * pscale / xstp[j]
        end
    end
    return inz
end

end # module 