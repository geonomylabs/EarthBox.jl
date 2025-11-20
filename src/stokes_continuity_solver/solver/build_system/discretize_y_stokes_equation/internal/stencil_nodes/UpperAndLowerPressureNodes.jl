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

    # Coefficient for the upper node p'U
    iprU = ivy + 1
    Lv[inz] = pscale/ystpc[i+1]
    Li[inz] = ivy
    Lj[inz] = iprU
    Lii[inz] = i
    Ljj[inz] = j
    inz = inz + 1

    # Bottom P node
    iprD = ivy + 4
    Lv[inz] = -pscale/ystpc[i+1]
    Li[inz] = ivy
    Lj[inz] = iprD
    Lii[inz] = i
    Ljj[inz] = j
    inz = inz + 1

    return inz
end

end # module 