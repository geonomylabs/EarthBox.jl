module UpdateStaggeredGrids

import EarthBox.ModelDataContainer.Grids2dContainer: Grids
import EarthBox.ModelDataContainer.Grids3dContainer: Grids3d
import ..PressureGridTools
import ..VxGridTools
import ..VyGridTools
import ..VzGridTools
import ..SxyGridTools
import ..SxzGridTools
import ..SyzGridTools

""" Update 2D staggered vx, vy and pressure grids and grid steps

Staggered Grid

    vx       vx       vx

vy  T---vy---T---vy---T   vy
    |        |        |
    vx   P   vx   P   vx
    |        |        |
vy  T---vy---T---vy---T   vy
    |        |        |
    vx   P   vx   P   vx
    |        |        |
vy  T---vy---T---vy---T   vy

    vx       vx       vx
"""
function update_staggered_grids!(grids::Grids)::Nothing
    # Calculate coordinates before spacing
    VyGridTools.calculate_xgrid_coordinates!(grids)
    VxGridTools.calculate_ygrid_coordinates!(grids)
    # Use coordinates to calculate spacing
    VyGridTools.calculate_xspacing!(grids)
    VxGridTools.calculate_yspacing!(grids)
    # For pressure calculate spacing using the basic grid first
    PressureGridTools.calculate_grid_spacing!(grids)
    # Calculate coordinates using spacing
    PressureGridTools.calculate_grid_coordinates!(grids)
    return nothing
end

""" Update 3D staggered vx, vy and vz grids and grid steps
"""
function update_staggered_grids!(grids3d::Grids3d)::Nothing
    # Calculate coordinates before spacing
    VyGridTools.calculate_ygrid_coordinates!(grids3d)
    VyGridTools.calculate_xgrid_coordinates!(grids3d)
    VyGridTools.calculate_zgrid_coordinates!(grids3d)

    VxGridTools.calculate_ygrid_coordinates!(grids3d)
    VxGridTools.calculate_xgrid_coordinates!(grids3d)
    VxGridTools.calculate_zgrid_coordinates!(grids3d)

    VzGridTools.calculate_xgrid_coordinates!(grids3d)
    VzGridTools.calculate_ygrid_coordinates!(grids3d)
    VzGridTools.calculate_zgrid_coordinates!(grids3d)

    SxyGridTools.calculate_xgrid_coordinates!(grids3d)
    SxyGridTools.calculate_ygrid_coordinates!(grids3d)
    SxyGridTools.calculate_zgrid_coordinates!(grids3d)

    SxzGridTools.calculate_xgrid_coordinates!(grids3d)
    SxzGridTools.calculate_ygrid_coordinates!(grids3d)
    SxzGridTools.calculate_zgrid_coordinates!(grids3d)

    SyzGridTools.calculate_xgrid_coordinates!(grids3d)
    SyzGridTools.calculate_ygrid_coordinates!(grids3d)
    SyzGridTools.calculate_zgrid_coordinates!(grids3d)

    # Use coordinates to calculate spacing
    VyGridTools.calculate_yspacing!(grids3d)
    VyGridTools.calculate_xspacing!(grids3d)
    VyGridTools.calculate_zspacing!(grids3d)

    VxGridTools.calculate_yspacing!(grids3d)
    VxGridTools.calculate_xspacing!(grids3d)
    VxGridTools.calculate_zspacing!(grids3d)

    VzGridTools.calculate_yspacing!(grids3d)
    VzGridTools.calculate_xspacing!(grids3d)
    VzGridTools.calculate_zspacing!(grids3d)

    SxyGridTools.calculate_xspacing!(grids3d)
    SxyGridTools.calculate_yspacing!(grids3d)
    SxyGridTools.calculate_zspacing!(grids3d)

    SxzGridTools.calculate_xspacing!(grids3d)
    SxzGridTools.calculate_yspacing!(grids3d)
    SxzGridTools.calculate_zspacing!(grids3d)   

    SyzGridTools.calculate_xspacing!(grids3d)
    SyzGridTools.calculate_yspacing!(grids3d)
    SyzGridTools.calculate_zspacing!(grids3d)

    # For pressure calculate spacing using the basic grid first
    PressureGridTools.calculate_grid_spacing!(grids3d)
    # Calculate coordinates using spacing
    PressureGridTools.calculate_grid_coordinates!(grids3d)
    return nothing
end

end # module 