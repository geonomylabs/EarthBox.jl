module BottomVelocityYNode

import EarthBox.Domain: basic_cell_in_last_two_bottom_rows
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
    ynum = build_data.grid.ynum
    ystp = build_data.grid.ystp
    ystpc = build_data.grid.ystpc
    etan = build_data.etan
    bbottomy = build_data.bc.bbottomy
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    if !basic_cell_in_last_two_bottom_rows(i, ynum)
        ivyD = ivy + 3
        Lv[inz] = 2.0 * etan[i+1, j] / ystp[i+1] / ystpc[i+1]
        Li[inz] = ivy
        Lj[inz] = ivyD
        Lii[inz] = i
        Ljj[inz] = j
        inz = inz + 1
    else
        # Apply boundary conditions since vyD is a ghost boundary node.
        # This leads to new terms for the coefficient for vyC and the rhs term.
        Lv[inz_c] = (
            Lv[inz_c] +
            bbottomy[j+1, 2] * 2.0 * etan[i+1, j] / ystp[i+1] / ystpc[i+1]
        )
        R[ivy] = (
            R[ivy] -
            bbottomy[j+1, 1] * 2.0 * etan[i+1, j] / ystp[i+1] / ystpc[i+1]
        )
    end
    return inz
end

end # module 