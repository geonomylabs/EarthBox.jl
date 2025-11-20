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

    # Left P node
    iprL = ivx + 2
    Lv[inz] = pscale / xstpc[j+1]
    Li[inz] = ivx
    Lj[inz] = iprL
    Lii[inz] = i
    Ljj[inz] = j
    #if inz < 1000
    #    println("LeftAndRightPressureNodes P1: inz = $inz, Lv = $(Lv[inz])")
    #end
    
    inz = inz + 1


    # Right P node
    iprR = ivx + 2 + hshift_to_vxR
    Lv[inz] = -pscale / xstpc[j+1]
    Li[inz] = ivx
    Lj[inz] = iprR
    Lii[inz] = i
    Ljj[inz] = j
    #if inz < 1000
    #    println("LeftAndRightPressureNodes P2: inz = $inz, Lv = $(Lv[inz])")
    #end

    inz = inz + 1


    return inz
end

end # module 