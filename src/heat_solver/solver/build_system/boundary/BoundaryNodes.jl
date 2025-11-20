module BoundaryNodes

import EarthBox.Domain: basic_node_on_upper_inside_boundary
import EarthBox.Domain: basic_node_on_lower_inside_boundary
import EarthBox.Domain: basic_node_on_left_boundary
import EarthBox.Domain: basic_node_on_right_boundary
import ..BuildStructs: CellIndices
import ..BuildStructs: HeatBuildData

""" Apply heat equation stencil for boundary node (i,j).

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
    build_data::HeatBuildData
)::Int
    i = cell_indices.i
    j = cell_indices.j
    xnum = build_data.grid.xnum
    ynum = build_data.grid.ynum
    if basic_node_on_upper_inside_boundary(i, j, xnum)
        inz = coefficients_for_upper_boundary(inz, cell_indices, build_data)
    end
    if basic_node_on_lower_inside_boundary(i, j, xnum, ynum)
        inz = coefficients_for_lower_boundary(inz, cell_indices, build_data)
    end
    if basic_node_on_left_boundary(j)
        inz = coefficients_for_left_boundary(inz, cell_indices, build_data)
    end
    if basic_node_on_right_boundary(j, xnum)
        inz = coefficients_for_right_boundary(inz, cell_indices, build_data)
    end
    return inz
end

""" Coefficients for upper boundary nodes.

Upper boundary: tk[i,j] = btop(1) + btop(2)*tk(i+1,j)

# Arguments
- `inz::Int`: Number of non-zero matrix elements
- `cell_indices::CellIndices`: Cell indices
- `build_data::HeatBuildData`: Build data

# Returns
- `inz::Int`: Updated number of non-zero matrix elements
"""
function coefficients_for_upper_boundary(
    inz::Int,
    cell_indices::CellIndices,
    build_data::HeatBuildData
)::Int
    i = cell_indices.i
    j = cell_indices.j
    itk = cell_indices.itk
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R
    btop = build_data.bc.btop

    # Right part
    Rval = btop[j, 1]
    R[itk] = Rval

    Lval1 = 1.0
    Lv[inz] = Lval1
    Li[inz] = itk
    Lj[inz] = itk
    Lii[inz] = i
    Ljj[inz] = j
    inz += 1

    Lval2 = -btop[j, 2]
    Lv[inz] = Lval2
    Li[inz] = itk
    Lj[inz] = itk + 1
    Lii[inz] = i
    Ljj[inz] = j
    inz += 1

    return inz
end

""" Coefficients for lower boundary nodes.

Lower boundary: tk[i,j] = bbottom(2) + bbottom(3)*tk(i-1,j)

# Arguments
- `inz::Int`: Number of non-zero matrix elements
- `cell_indices::CellIndices`: Cell indices
- `build_data::HeatBuildData`: Build data

# Returns
- `inz::Int`: Updated number of non-zero matrix elements
"""
function coefficients_for_lower_boundary(
    inz::Int,
    cell_indices::CellIndices,
    build_data::HeatBuildData
)::Int
    i = cell_indices.i
    j = cell_indices.j
    itk = cell_indices.itk
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R
    bbottom = build_data.bc.bbottom

    # Right part
    Rval = bbottom[j, 1]
    R[itk] = Rval

    # Left part: 1*tk[i,j] - bbottom(2)*tk(i-1,j)
    Lval1 = 1.0
    Lv[inz] = Lval1
    Li[inz] = itk
    Lj[inz] = itk
    Lii[inz] = i
    Ljj[inz] = j
    inz += 1

    Lval2 = -bbottom[j, 2]
    Lv[inz] = Lval2
    Li[inz] = itk
    Lj[inz] = itk - 1
    Lii[inz] = i
    Ljj[inz] = j
    inz += 1

    return inz
end

""" Coefficients for left boundary nodes.

Left boundary: tk[i,j] = bleft(1) + bleft(2)*tk[i,j+1]

# Arguments
- `inz::Int`: Number of non-zero matrix elements
- `cell_indices::CellIndices`: Cell indices
- `build_data::HeatBuildData`: Build data

# Returns
- `inz::Int`: Updated number of non-zero matrix elements
"""
function coefficients_for_left_boundary(
    inz::Int,
    cell_indices::CellIndices,
    build_data::HeatBuildData
)::Int
    i = cell_indices.i
    j = cell_indices.j
    itk = cell_indices.itk
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R
    bleft = build_data.bc.bleft
    ynum = build_data.grid.ynum

    # Right part
    Rval = bleft[i, 1]
    R[itk] = Rval

    # Left part: 1*tk[i,j] - bleft(2)*tk[i,j+1] = 0
    Lval1 = 1.0
    Lv[inz] = Lval1
    Li[inz] = itk
    Lj[inz] = itk
    Lii[inz] = i
    Ljj[inz] = j
    inz += 1

    Lval2 = -bleft[i, 2]
    Lv[inz] = Lval2
    Li[inz] = itk
    Lj[inz] = itk + ynum
    Lii[inz] = i
    Ljj[inz] = j
    inz += 1

    return inz
end

""" Coefficients for right boundary nodes.

Right boundary: tk[i,j] = bright(1) + bright(2)*tk[i,j-1]

# Arguments
- `inz::Int`: Number of non-zero matrix elements
- `cell_indices::CellIndices`: Cell indices
- `build_data::HeatBuildData`: Build data

# Returns
- `inz::Int`: Updated number of non-zero matrix elements
"""
function coefficients_for_right_boundary(
    inz::Int,
    cell_indices::CellIndices,
    build_data::HeatBuildData
)::Int
    i = cell_indices.i
    j = cell_indices.j
    itk = cell_indices.itk
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R
    bright = build_data.bc.bright
    ynum = build_data.grid.ynum

    # Right part
    Rval = bright[i, 1]
    R[itk] = Rval

    # Left part: 1*tk[i,j]-bright(2)*tk[i,j-1]
    Lval1 = 1.0
    Lv[inz] = Lval1
    Li[inz] = itk
    Lj[inz] = itk
    Lii[inz] = i
    Ljj[inz] = j
    inz += 1

    Lval2 = -bright[i, 2]
    Lv[inz] = Lval2
    Li[inz] = itk
    Lj[inz] = itk - ynum
    Lii[inz] = i
    Ljj[inz] = j
    inz += 1

    return inz
end

end # module 