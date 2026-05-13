module FrontAndBackPressureNodes

import ....BuildStructs: StokesBuildData
import ....BuildStructs: CellIndices

function coefficients_and_rhs_term(
    inz::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Int
    i = cell_indices.i
    j = cell_indices.j
    k = cell_indices.k
    ivz = cell_indices.ivz
    zstpc = build_data.grid.zstpc
    dshift_to_vxF = build_data.grid.dshift_to_vxF
    pscale = build_data.pscale
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv

    # Front P node — pressure of cell (i, j, k). Within-cell offset to p'
    # from vz is +1 in 3D.
    iprF = ivz + 1
    Lv[inz] = pscale / zstpc[k+1]
    Li[inz] = ivz
    Lj[inz] = iprF
    Lii[inz] = i
    Ljj[inz] = j
    inz = inz + 1

    # Back P node — pressure of cell (i, j, k+1).
    iprB = ivz + 1 + dshift_to_vxF
    Lv[inz] = -pscale / zstpc[k+1]
    Li[inz] = ivz
    Lj[inz] = iprB
    Lii[inz] = i
    Ljj[inz] = j
    inz = inz + 1

    return inz
end

end # module
