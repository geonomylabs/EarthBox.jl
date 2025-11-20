module BoundaryCells

import EarthBox.Domain: basic_cell_along_left_boundary
import EarthBox.Domain: basic_cell_along_right_boundary
import EarthBox.Domain: basic_cell_along_upper_boundary
import ...BuildStructs: StokesBuildData
import ...BuildStructs: CellIndices

""" Calculate coefficients and rhs for p unknowns in boundary cells.

Internal vs boundary cells depend on pressure boundary condition type
(pfirst_mode). If pfirst_mode = 0, then the upper left corner cell is the
boundary cell while all others are internal cells for pressure. If
pfirst_mode = 1 then the top and bottom rows of cells are boundary cells
while all others are internal. If pfirst_mode = 2 the side columns of cells
are boundary cells while all others are internal.

The approach described in Figure 4 of the main module docstring is used to
discretize the continuity equation.

# Arguments
- `inz::Int`: 
    - Index of non-zero
      matrix arrays (Lii, Ljj, Li, Lj and Lv).

- `cell_indices::CellIndices`: 
    - Cell index information for the current cell.

- `build_data::StokesBuildData`: 
    - Build data.

# Returns
- `inz::Int`: Index of non-zero
  matrix arrays (Lii, Ljj, Li, Lj and Lv).
"""
function coefficients_and_rhs_terms(
    inz::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Int
    if build_data.bc.pressure_bc_mode == 0
        inz = pressure_defined_in_upper_left_corner(inz, cell_indices, build_data)
    elseif build_data.bc.pressure_bc_mode == 1
        inz = pressure_defined_at_top_and_bottom(inz, cell_indices, build_data)
    elseif build_data.bc.pressure_bc_mode == 2
        inz = pressure_defined_at_left_and_right(inz, cell_indices, build_data)
    end
    return inz
end

""" Calculate coefficients and rhs for p unknowns in corner cell.

This function is used for cases where the pressure boundary condition is
defined in the upper left cell of the model domain (pressure_bc_mode = 0).

The approach described in Figure 4 of the main module docstring is used to
discretize the continuity equation.

# Arguments
- `inz::Int`: Current index of non-zero matrix arrays.
- `cell_indices::CellIndices`: Cell index information for the current cell.
- `build_data::StokesBuildData`: Build data.

# Returns
- `inz::Int`: Updated index of non-zero matrix arrays.
"""
function pressure_defined_in_upper_left_corner(
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
    pscale = build_data.pscale
    pressure_bc = build_data.bc.pressure_bc

    # value is 1 scaled up for uniformity
    Lv[inz] = 2.0 * pscale / (xstpavg + ystpavg)
    Li[inz] = ipr
    Lj[inz] = ipr
    Lii[inz] = i
    Ljj[inz] = j
    inz = inz + 1
    # Since we are solving for scaled pressure. The pressure
    # boundary condition is prbc/pscale leading to the
    # following expression without pscale
    R[ipr] = 2.0 * pressure_bc / (xstpavg + ystpavg)
    return inz
end

""" Calculate coefficients and rhs for p unknowns in corner cell.

This function is used for cases where the pressure boundary condition is
defined at the top and bottom of the model domain (pressure_bc_mode = 1).
This option is used for the channel flow benchmark.

The approach described in Figure 4 of the main module docstring is used to
discretize the continuity equation.

# Arguments
- `inz::Int`: Current index of non-zero matrix arrays.
- `cell_indices::CellIndices`: Cell index information for the current cell.
- `build_data::StokesBuildData`: Build data.

# Returns
- `inz::Int`: Updated index of non-zero matrix arrays.
"""
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
    pscale = build_data.pscale
    pressure_bc = build_data.bc.pressure_bc

    Lv[inz] = 2.0 * pscale / (xstpavg + ystpavg)
    Li[inz] = ipr
    Lj[inz] = ipr
    Lii[inz] = i
    Ljj[inz] = j
    inz = inz + 1
    if basic_cell_along_upper_boundary(i)
        # Since we are solving for scaled pressure. The
        # pressure boundary condition is prbc/pscale leading
        # to the following expression without pscale
        R[ipr] = 2.0 * pressure_bc / (xstpavg + ystpavg)
    else
        R[ipr] = 0.0
    end
    return inz
end

""" Calculate coefficients and rhs for p unknowns in corner cell.

This function is used for cases where the pressure boundary condition is
defined at the left and right of the model domain (pressure_bc_mode = 2).

The approach described in Figure 4 of the main module docstring is used to
discretize the continuity equation.

# Arguments
- `inz::Int`: Current index of non-zero matrix arrays.
- `cell_indices::CellIndices`: Cell index information for the current cell.
- `build_data::StokesBuildData`: Build data.

# Returns
- `inz::Int`: Updated index of non-zero matrix arrays.
"""
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
    pscale = build_data.pscale
    pressure_bc = build_data.bc.pressure_bc

    Lv[inz] = 2.0 * pscale / (xstpavg + ystpavg)
    Li[inz] = ipr
    Lj[inz] = ipr
    Lii[inz] = i
    Ljj[inz] = j
    inz = inz + 1
    if basic_cell_along_left_boundary(j)
        # Since we are solving for scaled pressure. The
        # pressure boundary condition is prbc/pscale leading
        # to the following expression without pscale
        R[ipr] = 2.0 * pressure_bc / (xstpavg + ystpavg)
    else
        R[ipr] = 0.0
    end
    return inz
end

end # module 