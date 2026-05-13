module UpperAndLowerPressureNodes

import ....BuildStructs: StokesBuildData
import ....BuildStructs: CellIndices

function coefficients_and_rhs_term(
    inz::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Int
    i = cell_indices.i
    j = cell_indices.j
    ivy = cell_indices.ivy
    ystpc = build_data.grid.ystpc
    pscale = build_data.pscale
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv

    # Upper P node — pressure of cell (i, j, k). Within-cell offset to p'
    # is +2 from ivy in 3D.
    iprU = ivy + 2
    Lv[inz] = pscale / ystpc[i+1]
    Li[inz] = ivy
    Lj[inz] = iprU
    Lii[inz] = i
    Ljj[inz] = j
    inz = inz + 1

    # Lower P node — pressure of cell (i+1, j, k).
    iprD = ivy + 6
    Lv[inz] = -pscale / ystpc[i+1]
    Li[inz] = ivy
    Lj[inz] = iprD
    Lii[inz] = i
    Ljj[inz] = j
    inz = inz + 1

    return inz
end

end # module
