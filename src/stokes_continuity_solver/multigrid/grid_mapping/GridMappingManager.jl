module GridMappingManager

import EarthBox.Arrays.ArrayTypes.ScalarArray3D: ScalarArray3DState
import EarthBox.Arrays.ArrayTypes.ScalarArray2D: ScalarArray2DState
import EarthBox.ModelDataContainer.Grids3dContainer: Grids3d
import EarthBox.ModelDataContainer.Grids2dContainer: Grids
import EarthBox.Interpolation.MarkerGridMapping: get_index_of_left_node

const DEBUG = false

struct GridMapping
    IULBy::Array{Int64,3}
    JULBx::Array{Int64,3}
    KULBz::Array{Int64,3}
    DyULB::Array{Float64,3}
    DxULB::Array{Float64,3}
    DzULB::Array{Float64,3}
end

function GridMapping(
    scalar_array3d::ScalarArray3DState
)::GridMapping
    # Indices (i,j,) of upper-left (UL) node of cell in grid b that contains grid a node
    IULBy = zeros(Int64, size(scalar_array3d.array))
    JULBx = zeros(Int64, size(scalar_array3d.array))
    KULBz = zeros(Int64, size(scalar_array3d.array))
    # Normalized distance between grid a node and upper-left (UL) node of grid b cell
    DyULB = zeros(Float64, size(scalar_array3d.array))
    DxULB = zeros(Float64, size(scalar_array3d.array))
    DzULB = zeros(Float64, size(scalar_array3d.array))
    return GridMapping(IULBy, JULBx, KULBz, DxULB, DyULB, DzULB)
end

struct GridMapping2d
    IULy::Array{Int64,2}
    JULx::Array{Int64,2}
    DyUL::Array{Float64,2}
    DxUL::Array{Float64,2}
end

function GridMapping2d(
    scalar_array2d::ScalarArray2DState
)::GridMapping2d
    # Indices (i,j,) of upper-left (UL) node of cell in grid b that contains grid a node
    IULy = zeros(Int64, size(scalar_array2d.array))
    JULx = zeros(Int64, size(scalar_array2d.array))
    # Normalized distance between grid a node and upper-left (UL) node of grid b cell
    DyUL = zeros(Float64, size(scalar_array2d.array))
    DxUL = zeros(Float64, size(scalar_array2d.array))
    return GridMapping2d(IULy, JULx, DyUL, DxUL)
end

struct MappingGroup
    vx_map::GridMapping
    vy_map::GridMapping
    vz_map::GridMapping
    pr_map::GridMapping
    etaxy_map::GridMapping
    etaxz_map::GridMapping
    etayz_map::GridMapping
end

function MappingGroup(
    vx::ScalarArray3DState,
    vy::ScalarArray3DState,
    vz::ScalarArray3DState,
    pr::ScalarArray3DState,
    etaxy::ScalarArray3DState,
    etaxz::ScalarArray3DState,
    etayz::ScalarArray3DState,
)::MappingGroup
    vx_map = GridMapping(vx)
    vy_map = GridMapping(vy)
    vz_map = GridMapping(vz)
    pr_map = GridMapping(pr)
    etaxy_map = GridMapping(etaxy)
    etaxz_map = GridMapping(etaxz)
    etayz_map = GridMapping(etayz)
    return MappingGroup(vx_map, vy_map, vz_map, pr_map, etaxy_map, etaxz_map, etayz_map)
end

struct MappingGroup2d
    vx_map::GridMapping2d
    vy_map::GridMapping2d
    pr_map::GridMapping2d
    etas_map::GridMapping2d
end

function MappingGroup2d(
    vx::ScalarArray2DState,
    vy::ScalarArray2DState,
    pr::ScalarArray2DState,
    etas::ScalarArray2DState,
)::MappingGroup2d
    vx_map = GridMapping2d(vx)
    vy_map = GridMapping2d(vy)
    pr_map = GridMapping2d(pr)
    etas_map = GridMapping2d(etas)
    return MappingGroup2d(vx_map, vy_map, pr_map, etas_map)
end

""" Calculate mapping group for mapping form grid a to grid b 
"""
function calculate_mapping_group!(
    mapping_group::MappingGroup,
    grid_a::Grids3d,
    grid_b::Grids3d,
)::Nothing
    calculate_grid_to_grid_mapping!(
        grid_a.arrays.staggered_vx.gridy_vx.array, grid_a.arrays.staggered_vx.gridx_vx.array, 
        grid_a.arrays.staggered_vx.gridz_vx.array, 
        grid_b.arrays.staggered_vx.gridy_vx.array, grid_b.arrays.staggered_vx.gridx_vx.array, 
        grid_b.arrays.staggered_vx.gridz_vx.array, 
        mapping_group.vx_map
        )
    calculate_grid_to_grid_mapping!(
        grid_a.arrays.staggered_vy.gridy_vy.array, grid_a.arrays.staggered_vy.gridx_vy.array, 
        grid_a.arrays.staggered_vy.gridz_vy.array, 
        grid_b.arrays.staggered_vy.gridy_vy.array, grid_b.arrays.staggered_vy.gridx_vy.array, 
        grid_b.arrays.staggered_vy.gridz_vy.array, 
        mapping_group.vy_map
        )
    calculate_grid_to_grid_mapping!(
        grid_a.arrays.staggered_vz.gridy_vz.array, grid_a.arrays.staggered_vz.gridx_vz.array, 
        grid_a.arrays.staggered_vz.gridz_vz.array, 
        grid_b.arrays.staggered_vz.gridy_vz.array, grid_b.arrays.staggered_vz.gridx_vz.array, 
        grid_b.arrays.staggered_vz.gridz_vz.array, 
        mapping_group.vz_map
        )
    calculate_grid_to_grid_mapping!(
        grid_a.arrays.pressure.gridy_pr.array, grid_a.arrays.pressure.gridx_pr.array, 
        grid_a.arrays.pressure.gridz_pr.array, 
        grid_b.arrays.pressure.gridy_pr.array, grid_b.arrays.pressure.gridx_pr.array, 
        grid_b.arrays.pressure.gridz_pr.array, 
        mapping_group.pr_map
        )
    calculate_grid_to_grid_mapping!(
        grid_a.arrays.staggered_sxy.gridy_sxy.array, grid_a.arrays.staggered_sxy.gridx_sxy.array, 
        grid_a.arrays.staggered_sxy.gridz_sxy.array, 
        grid_b.arrays.staggered_sxy.gridy_sxy.array, grid_b.arrays.staggered_sxy.gridx_sxy.array, 
        grid_b.arrays.staggered_sxy.gridz_sxy.array, 
        mapping_group.etaxy_map
        )
    calculate_grid_to_grid_mapping!(
        grid_a.arrays.staggered_sxz.gridy_sxz.array, grid_a.arrays.staggered_sxz.gridx_sxz.array, 
        grid_a.arrays.staggered_sxz.gridz_sxz.array, 
        grid_b.arrays.staggered_sxz.gridy_sxz.array, grid_b.arrays.staggered_sxz.gridx_sxz.array, 
        grid_b.arrays.staggered_sxz.gridz_sxz.array, 
        mapping_group.etaxz_map
        )
    calculate_grid_to_grid_mapping!(
        grid_a.arrays.staggered_syz.gridy_syz.array, grid_a.arrays.staggered_syz.gridx_syz.array, 
        grid_a.arrays.staggered_syz.gridz_syz.array, 
        grid_b.arrays.staggered_syz.gridy_syz.array, grid_b.arrays.staggered_syz.gridx_syz.array, 
        grid_b.arrays.staggered_syz.gridz_syz.array, 
        mapping_group.etayz_map
        )
end

"""  Calculate mapping info on grid a for mapping from from grid a to grid b.

Mapping is calculated relative to the upper-left-back (ULB) node of the grid b cell that contains 
the grid a node. Upper is in the negative y-direction where y is positive in the downward direction.
Left refers to the negative x-direction. Right refers to the positive x-direction. Back refers to the 
negative z-direction. Front refers to the positive z-direction.
    
           ___________
          /|         /|         UP
         / |        / |
 ULB--> /__|_______/__| Front   Down 
        |          |  /
        |          | /
        |__________|/ Back
     Left          Right        

"""
function calculate_grid_to_grid_mapping!(
    gridy_a::Array{Float64,1},
    gridx_a::Array{Float64,1},
    gridz_a::Array{Float64,1},
    gridy_b::Array{Float64,1},
    gridx_b::Array{Float64,1},
    gridz_b::Array{Float64,1},
    grid_mapping::GridMapping
)::Nothing
    # Check if grid sizes match and throw error if they are the same
    if size(gridy_a, 1) == size(gridy_b, 1) || 
       size(gridx_a, 1) == size(gridx_b, 1) ||
       size(gridz_a, 1) == size(gridz_b, 1)
        error("Grid sizes should not be identical for mapping between grids")
    end

    # unpack grid mapping arrays
    IULBy = grid_mapping.IULBy
    JULBx = grid_mapping.JULBx
    KULBz = grid_mapping.KULBz
    DyULB = grid_mapping.DyULB
    DxULB = grid_mapping.DxULB
    DzULB = grid_mapping.DzULB

    ynuma = size(gridy_a, 1)
    xnuma = size(gridx_a, 1)
    znuma = size(gridz_a, 1)
    ynumb = size(gridy_b, 1)
    xnumb = size(gridx_b, 1)
    znumb = size(gridz_b, 1)

    if DEBUG
        println("xnuma = $(xnuma), xnumb = $(xnumb)")
        println("gridx_a = $(gridx_a)")
        println("gridx_b = $(gridx_b)")
        println("znuma = $(znuma), znumb = $(znumb)")
        println("gridz_a = $(gridz_a)")
        println("gridz_b = $(gridz_b)")
        println("ynuma = $(ynuma), ynumb = $(ynumb)")
        println("gridy_a = $(gridy_a)")
        println("gridy_b = $(gridy_b)")
    end

    for k in 1:znuma
        for j in 1:xnuma
            for i in 1:ynuma
                # Coordinates of node on grid a
                ya = gridy_a[i]
                xa = gridx_a[j]
                za = gridz_a[k]
                # Calculate indices (iULB, jULB, kULB) of upper-left-back (ULB) node of cell in 
                # grid b that contains grid a node
                iULBy = get_index_of_left_node(ya, gridy_b, ynumb)
                jULBx = get_index_of_left_node(xa, gridx_b, xnumb)
                kULBz = get_index_of_left_node(za, gridz_b, znumb)
                IULBy[i,j,k] = iULBy
                JULBx[i,j,k] = jULBx
                KULBz[i,j,k] = kULBz
                # Coordinates of upper-left-back (ULB) node of cell in grid b
                ybULB = gridy_b[iULBy]
                xbULB = gridx_b[jULBx]
                zbULB = gridz_b[kULBz]
                # Height of grid b cell that contains the grid a node
                cell_height_b = gridy_b[iULBy+1] - ybULB
                # Width of grid b cell that contains the grid a node
                cell_width_b = gridx_b[jULBx+1] - xbULB
                # Depth of grid b cell that contains the grid a node
                cell_depth_b = gridz_b[kULBz+1] - zbULB
                # Normalized distance between grid a node and upper-left (UL) node of grid b cell
                dyULB = (ya - ybULB)/cell_height_b
                dxULB = (xa - xbULB)/cell_width_b
                dzULB = (za - zbULB)/cell_depth_b
                DxULB[i,j,k] = dxULB
                DyULB[i,j,k] = dyULB
                DzULB[i,j,k] = dzULB
            end
        end
    end
    
end

""" Calculate mapping group for mapping form grid a to grid b for 2d
"""
function calculate_mapping_group!(
    mapping_group::MappingGroup2d,
    grid_a::Grids,
    grid_b::Grids,
)::Nothing
    calculate_grid_to_grid_mapping!(
        grid_a.arrays.staggered_vx.gridy_vx.array, grid_a.arrays.basic.gridx_b.array, 
        grid_b.arrays.staggered_vx.gridy_vx.array, grid_b.arrays.basic.gridx_b.array, 
        mapping_group.vx_map
        )
    calculate_grid_to_grid_mapping!(
        grid_a.arrays.basic.gridy_b.array, grid_a.arrays.staggered_vy.gridx_vy.array, 
        grid_b.arrays.basic.gridy_b.array, grid_b.arrays.staggered_vy.gridx_vy.array, 
        mapping_group.vy_map
        )
    calculate_grid_to_grid_mapping!(
        grid_a.arrays.pressure.gridy_pr.array, grid_a.arrays.pressure.gridx_pr.array, 
        grid_b.arrays.pressure.gridy_pr.array, grid_b.arrays.pressure.gridx_pr.array, 
        mapping_group.pr_map
        )
    calculate_grid_to_grid_mapping!(
        grid_a.arrays.basic.gridy_b.array, grid_a.arrays.basic.gridx_b.array, 
        grid_b.arrays.basic.gridy_b.array, grid_b.arrays.basic.gridx_b.array, 
        mapping_group.etas_map
        )
end

"""  Calculate mapping info on grid a for mapping from from grid a to grid b fro 2d
"""
function calculate_grid_to_grid_mapping!(
    gridy_a::Array{Float64,1},
    gridx_a::Array{Float64,1},
    gridy_b::Array{Float64,1},
    gridx_b::Array{Float64,1},
    grid_mapping::GridMapping2d
)::Nothing
    # Check if grid sizes match and throw error if they are the same
    if size(gridy_a, 1) == size(gridy_b, 1) || 
       size(gridx_a, 1) == size(gridx_b, 1)
        error("Grid sizes should not be identical for mapping between grids")
    end

    # unpack grid mapping arrays
    IULy = grid_mapping.IULy
    JULx = grid_mapping.JULx
    DyUL = grid_mapping.DyUL
    DxUL = grid_mapping.DxUL

    ynuma = size(gridy_a, 1)
    xnuma = size(gridx_a, 1)
    ynumb = size(gridy_b, 1)
    xnumb = size(gridx_b, 1)

    for j in 1:xnuma
        for i in 1:ynuma
            # Coordinates of node on grid a
            ya = gridy_a[i]
            xa = gridx_a[j]
            # Calculate indices (iULB, jULB) of upper-left (UL) node of cell in 
            # grid b that contains grid a node
            iULy = get_index_of_left_node(ya, gridy_b, ynumb)
            jULx = get_index_of_left_node(xa, gridx_b, xnumb)
            IULy[i,j] = iULy
            JULx[i,j] = jULx
            # Coordinates of upper-left (UL) node of cell in grid b
            ybUL = gridy_b[iULy]
            xbUL = gridx_b[jULx]
            # Height of grid b cell that contains the grid a node
            cell_height_b = gridy_b[iULy+1] - ybUL
            # Width of grid b cell that contains the grid a node
            cell_width_b = gridx_b[jULx+1] - xbUL
            # Normalized distance between grid a node and upper-left (UL) node of grid b cell
            dyUL = (ya - ybUL)/cell_height_b
            dxUL = (xa - xbUL)/cell_width_b
            DxUL[i,j] = dxUL
            DyUL[i,j] = dyUL
        end
    end
end

end