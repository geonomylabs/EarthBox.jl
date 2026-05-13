module RightVelocityXNode

import EarthBox.Domain: basic_cell_in_last_two_rightmost_columns
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
    xnum = build_data.grid.xnum
    xstp = build_data.grid.xstp
    xstpc = build_data.grid.xstpc
    etan = build_data.etan
    brightx = build_data.bc.brightx
    hshift_to_vxR = build_data.grid.hshift_to_vxR
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    if !basic_cell_in_last_two_rightmost_columns(j, xnum)
        ivxR = ivx + hshift_to_vxR
        Lv[inz] = 2.0 * etan[i, j+1, k] / xstp[j+1] / xstpc[j+1]
        Li[inz] = ivx
        Lj[inz] = ivxR
        Lii[inz] = i
        Ljj[inz] = j
        inz = inz + 1
    else
        # vxR is the prescribed right-boundary velocity. Apply the
        # boundary condition equation and fold the result into vxC and
        # the rhs term.
        Lv[inz_c] = (
            Lv[inz_c] +
            brightx[i+1, k+1, 2] * 2.0 * etan[i, j+1, k] / xstp[j+1] / xstpc[j+1]
        )
        R[ivx] = (
            R[ivx] -
            brightx[i+1, k+1, 1] * 2.0 * etan[i, j+1, k] / xstp[j+1] / xstpc[j+1]
        )
    end
    return inz
end

end # module
