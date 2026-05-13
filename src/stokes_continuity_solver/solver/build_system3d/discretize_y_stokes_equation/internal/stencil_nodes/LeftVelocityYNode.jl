module LeftVelocityYNode

import EarthBox.Domain: basic_cell_along_left_boundary
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
    xstp = build_data.grid.xstp
    xstpc = build_data.grid.xstpc
    etas_xy = build_data.etas_xy
    blefty = build_data.bc.blefty
    hshift_to_vxR = build_data.grid.hshift_to_vxR
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    if !basic_cell_along_left_boundary(j)
        ivyL = ivy - hshift_to_vxR
        Lv[inz] = etas_xy[i+1, j, k] / xstpc[j] / xstp[j]
        Li[inz] = ivy
        Lj[inz] = ivyL
        Lii[inz] = i
        Ljj[inz] = j
        inz = inz + 1
    else
        # vyL is a left ghost boundary node.
        Lv[inz_c] = (
            Lv[inz_c] +
            blefty[i+1, k+1, 2] * etas_xy[i+1, j, k] / xstpc[j] / xstp[j]
        )
        R[ivy] = (
            R[ivy] -
            blefty[i+1, k+1, 1] * etas_xy[i+1, j, k] / xstpc[j] / xstp[j]
        )
    end
    return inz
end

end # module
