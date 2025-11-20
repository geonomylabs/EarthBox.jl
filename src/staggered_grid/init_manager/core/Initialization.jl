module Initialization

import EarthBox.ModelDataContainer.Grids2dContainer: Grids
import EarthBox.ModelDataContainer.Grids3dContainer: Grids3d
import ...TypesManager: Uniform
import ...TypesManager: SidesWithHighres
import ...TypesManager: SidesWithHighresAndSmoothing
import ...TypesManager: Ttype
import ...UpdateBasicAndStaggeredGrids

function initialize_grid_uniform!(grids::Grids)
    Uniform.Initialization.initialize!(grids)
    UpdateBasicAndStaggeredGrids.update!(grids)
end

function initialize_grid_uniform!(grids3d::Grids3d)
    Uniform.Initialization.initialize!(grids3d)
    UpdateBasicAndStaggeredGrids.update!(grids3d)
end

function initialize_grid_51x51_uniform_side_refinement!(grids::Grids)
    SidesWithHighres.Initialization.initialize!(grids)
    UpdateBasicAndStaggeredGrids.update!(grids)
end

function initialize_grid_51x51_smoothed_side_refinement!(grids::Grids)
    SidesWithHighresAndSmoothing.Initialization.initialize!(grids)
    UpdateBasicAndStaggeredGrids.update!(grids)
end

function initialize_grid_t_type_refined!(grids::Grids)
    Ttype.Initialization.initialize!(grids)
    UpdateBasicAndStaggeredGrids.update!(grids)
end

end # module Initialization

