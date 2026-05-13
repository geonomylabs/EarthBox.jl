module BoundaryCells

import EarthBox.Domain: basic_cell_along_left_boundary
import EarthBox.Domain: basic_cell_along_right_boundary
import EarthBox.Domain: basic_cell_along_upper_boundary
import ...Predicates3D: basic_cell_along_front_boundary
import ...BuildStructs: StokesBuildData
import ...BuildStructs: CellIndices

""" Calculate coefficients and rhs for p unknowns in 3D boundary cells.

Internal vs boundary cells depend on the pressure boundary-condition type
`pressure_bc_mode`:

  * mode = 0: only the upper-left-front corner cell (i = j = k = 1) is a
              boundary cell.
  * mode = 1: top and bottom y-face cells are boundary cells.
  * mode = 2: left and right x-face cells are boundary cells.
  * mode = 3: front and back z-face cells are boundary cells (new in 3D).

# Arguments
- `inz::Int`:
    - Index of non-zero matrix arrays (Lii, Ljj, Li, Lj and Lv).

- `cell_indices::CellIndices`:
    - Cell index information for the current cell.

- `build_data::StokesBuildData`:
    - Build data.

# Returns
- `inz::Int`: Updated index of non-zero matrix arrays.
"""
function coefficients_and_rhs_terms(
    inz::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Int
    mode = build_data.bc.pressure_bc_mode
    if mode == 0
        inz = pressure_defined_in_upper_left_front_corner(inz, cell_indices, build_data)
    elseif mode == 1
        inz = pressure_defined_at_top_and_bottom(inz, cell_indices, build_data)
    elseif mode == 2
        inz = pressure_defined_at_left_and_right(inz, cell_indices, build_data)
    elseif mode == 3
        inz = pressure_defined_at_front_and_back(inz, cell_indices, build_data)
    end
    return inz
end

function pressure_defined_in_upper_left_front_corner(
    inz::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Int
    i = cell_indices.i
    j = cell_indices.j
    ipr = cell_indices.ipr
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R
    xstpavg = build_data.grid.xstpavg
    ystpavg = build_data.grid.ystpavg
    zstpavg = build_data.grid.zstpavg
    pscale = build_data.pscale
    pressure_bc = build_data.bc.pressure_bc

    Lv[inz] = 2.0 * pscale / (xstpavg + ystpavg + zstpavg)
    Li[inz] = ipr
    Lj[inz] = ipr
    Lii[inz] = i
    Ljj[inz] = j
    inz = inz + 1
    # Solve for scaled pressure: prbc/pscale leaves the following without
    # pscale.
    R[ipr] = 2.0 * pressure_bc / (xstpavg + ystpavg + zstpavg)
    return inz
end

function pressure_defined_at_top_and_bottom(
    inz::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Int
    i = cell_indices.i
    j = cell_indices.j
    ipr = cell_indices.ipr
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R
    xstpavg = build_data.grid.xstpavg
    ystpavg = build_data.grid.ystpavg
    zstpavg = build_data.grid.zstpavg
    pscale = build_data.pscale
    pressure_bc = build_data.bc.pressure_bc

    Lv[inz] = 2.0 * pscale / (xstpavg + ystpavg + zstpavg)
    Li[inz] = ipr
    Lj[inz] = ipr
    Lii[inz] = i
    Ljj[inz] = j
    inz = inz + 1
    if basic_cell_along_upper_boundary(i)
        R[ipr] = 2.0 * pressure_bc / (xstpavg + ystpavg + zstpavg)
    else
        R[ipr] = 0.0
    end
    return inz
end

function pressure_defined_at_left_and_right(
    inz::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Int
    i = cell_indices.i
    j = cell_indices.j
    ipr = cell_indices.ipr
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R
    xstpavg = build_data.grid.xstpavg
    ystpavg = build_data.grid.ystpavg
    zstpavg = build_data.grid.zstpavg
    pscale = build_data.pscale
    pressure_bc = build_data.bc.pressure_bc

    Lv[inz] = 2.0 * pscale / (xstpavg + ystpavg + zstpavg)
    Li[inz] = ipr
    Lj[inz] = ipr
    Lii[inz] = i
    Ljj[inz] = j
    inz = inz + 1
    if basic_cell_along_left_boundary(j)
        R[ipr] = 2.0 * pressure_bc / (xstpavg + ystpavg + zstpavg)
    else
        R[ipr] = 0.0
    end
    return inz
end

function pressure_defined_at_front_and_back(
    inz::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Int
    i = cell_indices.i
    j = cell_indices.j
    k = cell_indices.k
    ipr = cell_indices.ipr
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R
    xstpavg = build_data.grid.xstpavg
    ystpavg = build_data.grid.ystpavg
    zstpavg = build_data.grid.zstpavg
    pscale = build_data.pscale
    pressure_bc = build_data.bc.pressure_bc

    Lv[inz] = 2.0 * pscale / (xstpavg + ystpavg + zstpavg)
    Li[inz] = ipr
    Lj[inz] = ipr
    Lii[inz] = i
    Ljj[inz] = j
    inz = inz + 1
    if basic_cell_along_front_boundary(k)
        R[ipr] = 2.0 * pressure_bc / (xstpavg + ystpavg + zstpavg)
    else
        R[ipr] = 0.0
    end
    return inz
end

end # module
