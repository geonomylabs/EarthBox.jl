module RightVelocityYNode

import EarthBox.Domain: basic_cell_along_right_boundary
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
    xstp = build_data.grid.xstp
    xstpc = build_data.grid.xstpc
    etas = build_data.etas
    brighty = build_data.bc.brighty
    hshift_to_vxR = build_data.grid.hshift_to_vxR
    xnum = build_data.grid.xnum
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    if !basic_cell_along_right_boundary(j, xnum)
        ivyR = ivy + hshift_to_vxR
        Lv[inz] = etas[i+1,j+1]/xstpc[j+1]/xstp[j]
        Li[inz] = ivy
        Lj[inz] = ivyR
        Lii[inz] = i
        Ljj[inz] = j
        inz = inz + 1
    else
        # Apply boundary conditions since vyR is a ghost boundary
        # node. This leads to new terms for the coefficient for
        # vyC and the rhs term.
        Lv[inz_c] += brighty[i+1,2]*etas[i+1,j+1]/xstpc[j+1]/xstp[j]
        R[ivy] -= brighty[i+1,1]*etas[i+1,j+1]/xstpc[j+1]/xstp[j]
    end
    return inz
end

end # module 