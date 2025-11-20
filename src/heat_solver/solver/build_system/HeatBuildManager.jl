""" Calculate non-zero coefficients and rhs vector for heat equation.

These functions formulate a discretized system of equations using a conservative 
finite-difference approximation 2D grid (Figure 1). Only the non-zero coefficients 
of matrix L are computed by this function. The matrix indices of the non-zero 
elements, the basic grid indices Lii and Ljj associated with each matrix element 
and the right-and-side vector R are also computed.

              | xstpc(j-1) |
        | xstp(j-1) |  xstp(j)  |
        +-------tk(i-1,j)-------+---------
        |       kt(i-1,j)       |
        |           |           |
        |           |           |    ystp(i-1)-----
        |           |           |
        |           |           |
    tk(i,j-1)----tk(i,j)-----tk(i,j+1)------  ystpc(i-1)
    kt(i,j-1)    kt(i,j)     kt(i,j+1)
        |       rhocp(i,j)      |
        |           |           |    ystp(i)-------
        |           |           |
        |           |           |
        +-------tk(i+1,j)-------+---------
                kt(i+1,j)

    Figure 1: Heat equation stencil
"""
module HeatBuildManager

include("utils/BuildStructs.jl")
include("boundary/BoundaryNodes.jl")
include("internal/InternalNodes.jl")

import EarthBox.ModelDataContainer: ModelData
import EarthBox.Domain: basic_node_on_boundary
import EarthBox: BuildSysTools
import EarthBox.BuildSysTools: SystemVectors
import .BuildStructs: CellIndices
import .BuildStructs: HeatBuildData
import .BuildStructs: HeatStencilData
import .InternalNodes
import .BoundaryNodes

function build_system_of_equations(model::ModelData)::SystemVectors
    """ Calculate non-zero coefficients and rhs vector for heat equation.

    Input Container
    ---------------
    model: ModelData
        Model interface object used to store parameters and arrays.

    Returns
    -------
    inz: Int64
        Number of non-zero matrix elements.
    """
    build_data = HeatBuildData(model)
    xnum = build_data.grid.xnum
    ynum = build_data.grid.ynum
    inz = 1
    for j in 1:xnum
        for i in 1:ynum
            itk = (j-1)*ynum + i
            cell_indices = CellIndices(i, j, itk)
            inz = apply_stencils(inz, cell_indices, build_data)
        end
    end
    # Subtract 1 from inz to account for the initial value of 1.
    return BuildSysTools.clean_non_zero_arrays!(inz-1, build_data.system_vectors)
end 

function apply_stencils(
    inz::Int64, 
    cell_indices::CellIndices, 
    build_data::HeatBuildData
)::Int64
    xnum = build_data.grid.xnum
    ynum = build_data.grid.ynum
    if basic_node_on_boundary(cell_indices.i, cell_indices.j, xnum, ynum)
        inz = BoundaryNodes.apply_stencil(inz, cell_indices, build_data)
    else
        inz = InternalNodes.apply_stencil(inz, cell_indices, build_data)
    end
    return inz
end

end # module