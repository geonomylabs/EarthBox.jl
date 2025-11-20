module GridUpdate

import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer.Grids3dContainer: Grids3d
import EarthBox.ModelDataContainer.Grids2dContainer: Grids
import ....UpdateBasicAndStaggeredGrids
import ..Initialization: initialize!

""" Update t-type grid using current t-type grid parameters.
"""
function update_t_type_grid!(model::ModelData)
    check_t_type_refinement_parameters(model)
    initialize!(model.grids)
    UpdateBasicAndStaggeredGrids.update!(model.grids)
end

function update_t_type_grid!(grids::Grids)
    check_t_type_refinement_parameters(grids)
    initialize!(grids)
    UpdateBasicAndStaggeredGrids.update!(grids)
end

function check_t_type_refinement_parameters(model::ModelData)::Nothing
    check_t_type_refinement_parameters(model.grids)
    return nothing
end

function check_t_type_refinement_parameters(grids::Grids)::Nothing
    xnum       = grids.parameters.geometry.xnum.value
    ynum       = grids.parameters.geometry.ynum.value
    dx_highres = grids.parameters.refinement.dx_highres.value
    dy_highres = grids.parameters.refinement.dy_highres.value
    xo_highres = grids.parameters.refinement.xo_highres.value
    xf_highres = grids.parameters.refinement.xf_highres.value
    yf_highres = grids.parameters.refinement.yf_highres.value

    nnodes_highres_x = Int64(floor((xf_highres-xo_highres)/dx_highres)) + 1
    nnodes_highres_y = Int64(floor(yf_highres/dy_highres)) + 1

    buffer_x = 5
    buffer_y = 5
    min_nodes_x = nnodes_highres_x + buffer_x*2
    min_nodes_y = nnodes_highres_y + buffer_y
    if min_nodes_x  > xnum
        throw(ArgumentError(
            "T-type refinement: xnum is too small for high-resolution "
            *"refinement. Increase xnum to at least $min_nodes_x."
            )
        )
    end
    if nnodes_highres_y + buffer_y > ynum
        throw(ArgumentError(
            "T-type refinement: ynum is too small for the high-resolution "
            *"refinement. Increase ynum to at least $min_nodes_y."
            )
        )
    end

    return nothing
end

end

