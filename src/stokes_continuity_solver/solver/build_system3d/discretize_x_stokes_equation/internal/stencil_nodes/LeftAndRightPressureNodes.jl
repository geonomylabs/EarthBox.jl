module LeftAndRightPressureNodes

import ....BuildStructs: StokesBuildData
import ....BuildStructs: CellIndices

function coefficients_and_rhs_term(
    inz::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Int
    ivx = cell_indices.ivx
    i = cell_indices.i
    j = cell_indices.j
    xstpc = build_data.grid.xstpc
    hshift_to_vxR = build_data.grid.hshift_to_vxR
    pscale = build_data.pscale
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv

    # Left P node — pressure of cell (i, j, k). Within-cell offset to p' is
    # +3 in 3D (the 4th unknown).
    iprL = ivx + 3
    Lv[inz] = pscale / xstpc[j+1]
    Li[inz] = ivx
    Lj[inz] = iprL
    Lii[inz] = i
    Ljj[inz] = j
    inz = inz + 1

    # Right P node — pressure of cell (i, j+1, k).
    iprR = ivx + 3 + hshift_to_vxR
    Lv[inz] = -pscale / xstpc[j+1]
    Li[inz] = ivx
    Lj[inz] = iprR
    Lii[inz] = i
    Ljj[inz] = j
    inz = inz + 1

    return inz
end

end # module
