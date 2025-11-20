module InternalNodes

import ..BuildStructs: HeatStencilData
import ..BuildStructs: CellIndices
import ..BuildStructs: HeatBuildData

""" Apply heat equation stencil for internal node (i,j).

Temperature equation:
    RHO*CP*DT/Dt = -dqx/dx - dqy/dy + Ht

# Arguments
- `inz::Int`: Number of non-zero matrix elements
- `cell_indices::CellIndices`: Cell indices
- `build_data::HeatBuildData`: Build data

# Returns
- `inz::Int`: Updated number of non-zero matrix elements
"""
function apply_stencil(
    inz::Int,
    cell_indices::CellIndices,
    build_data::HeatBuildData,
)::Int
    stencil_data = HeatStencilData(build_data, cell_indices)
    update_right_part(cell_indices, build_data, stencil_data)
    inz = central_node(inz, cell_indices, build_data, stencil_data)
    inz = left_node(inz, cell_indices, build_data, stencil_data)
    inz = right_node(inz, cell_indices, build_data, stencil_data)
    inz = upper_node(inz, cell_indices, build_data, stencil_data)
    inz = lower_node(inz, cell_indices, build_data, stencil_data)
    return inz
end

function update_right_part(
    cell_indices::CellIndices,
    build_data::HeatBuildData,
    stencil_data::HeatStencilData,
)::Nothing
    i = cell_indices.i
    j = cell_indices.j
    itk = cell_indices.itk
    RT = build_data.rhs.RT
    R = build_data.rhs.R
    tkC = stencil_data.tkC
    rhocpC = stencil_data.rhocpC
    timestep = build_data.timestep

    Rval = RT[i, j] + tkC * rhocpC / timestep
    R[itk] = Rval
    return nothing
end

function central_node(
    inz::Int,
    cell_indices::CellIndices,
    build_data::HeatBuildData,
    stencil_data::HeatStencilData,
)::Int
    i = cell_indices.i
    j = cell_indices.j
    itk = cell_indices.itk
    ktL = stencil_data.ktL
    ktC = stencil_data.ktC
    ktR = stencil_data.ktR
    ktU = stencil_data.ktU
    ktD = stencil_data.ktD
    rhocpC = stencil_data.rhocpC
    dxL = stencil_data.dxL
    dxR = stencil_data.dxR
    dxC = stencil_data.dxC
    dyU = stencil_data.dyU
    dyD = stencil_data.dyD
    dyC = stencil_data.dyC
    timestep = build_data.timestep
    
    Lval1 = (
        rhocpC / timestep +
        ((ktL + ktC) / dxL + (ktC + ktR) / dxR) / 2.0 / dxC +
        ((ktU + ktC) / dyU + (ktC + ktD) / dyD) / 2.0 / dyC
    )

    build_data.system_vectors.Lv[inz] = Lval1
    build_data.system_vectors.Li[inz] = itk
    build_data.system_vectors.Lj[inz] = itk
    build_data.system_vectors.Lii[inz] = i
    build_data.system_vectors.Ljj[inz] = j

    inz += 1
    return inz
end

function left_node(
    inz::Int,
    cell_indices::CellIndices,
    build_data::HeatBuildData,
    stencil_data::HeatStencilData,
)::Int
    i = cell_indices.i
    j = cell_indices.j
    itk = cell_indices.itk
    ktC = stencil_data.ktC
    ktL = stencil_data.ktL
    dxL = stencil_data.dxL
    dxC = stencil_data.dxC
    ynum = build_data.grid.ynum

    Lval2 = -(ktL + ktC) / 2.0 / dxL / dxC

    build_data.system_vectors.Lv[inz] = Lval2
    build_data.system_vectors.Li[inz] = itk
    build_data.system_vectors.Lj[inz] = itk - ynum
    build_data.system_vectors.Lii[inz] = i
    build_data.system_vectors.Ljj[inz] = j

    inz += 1
    return inz
end

function right_node(
    inz::Int,
    cell_indices::CellIndices,
    build_data::HeatBuildData,
    stencil_data::HeatStencilData,
)::Int
    i = cell_indices.i
    j = cell_indices.j
    itk = cell_indices.itk
    ktC = stencil_data.ktC
    ktR = stencil_data.ktR
    dxR = stencil_data.dxR
    dxC = stencil_data.dxC
    ynum = build_data.grid.ynum

    Lval3 = -(ktC + ktR) / 2.0 / dxR / dxC

    build_data.system_vectors.Lv[inz] = Lval3
    build_data.system_vectors.Li[inz] = itk
    build_data.system_vectors.Lj[inz] = itk + ynum
    build_data.system_vectors.Lii[inz] = i
    build_data.system_vectors.Ljj[inz] = j

    inz += 1
    return inz
end

function upper_node(
    inz::Int,
    cell_indices::CellIndices,
    build_data::HeatBuildData,
    stencil_data::HeatStencilData,
)::Int
    i = cell_indices.i
    j = cell_indices.j
    itk = cell_indices.itk
    ktU = stencil_data.ktU
    ktC = stencil_data.ktC
    dyU = stencil_data.dyU
    dyC = stencil_data.dyC

    Lval4 = -(ktU + ktC) / 2.0 / dyU / dyC

    build_data.system_vectors.Lv[inz] = Lval4
    build_data.system_vectors.Li[inz] = itk
    build_data.system_vectors.Lj[inz] = itk - 1
    build_data.system_vectors.Lii[inz] = i
    build_data.system_vectors.Ljj[inz] = j

    inz += 1
    return inz
end

function lower_node(
    inz::Int,
    cell_indices::CellIndices,
    build_data::HeatBuildData,
    stencil_data::HeatStencilData,
)::Int
    i = cell_indices.i
    j = cell_indices.j
    itk = cell_indices.itk
    ktC = stencil_data.ktC
    ktD = stencil_data.ktD
    dyD = stencil_data.dyD
    dyC = stencil_data.dyC

    Lval5 = -(ktC + ktD) / 2.0 / dyD / dyC

    build_data.system_vectors.Lv[inz] = Lval5
    build_data.system_vectors.Li[inz] = itk
    build_data.system_vectors.Lj[inz] = itk + 1
    build_data.system_vectors.Lii[inz] = i
    build_data.system_vectors.Ljj[inz] = j

    inz += 1
    return inz
end

end # module 