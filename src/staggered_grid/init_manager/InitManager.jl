module InitManager

include("core/Initialization.jl")

import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer.Grids3dContainer: Grids3d
import EarthBox.ModelDataContainer.Grids2dContainer: Grids
import ..Options: option_names
import .Initialization
import ..TypesManager: Ttype

function initialize!(
    model::ModelData,
    ::Val{option_names.UniformGrid}
)
    initialize!(model.grids, Val(option_names.UniformGrid))
end

function initialize!(
    grids::Grids,
    ::Val{option_names.UniformGrid}
)
    Initialization.initialize_grid_uniform!(grids)
end

function initialize!(
    grids::Grids3d,
    ::Val{option_names.UniformGrid}
)
    Initialization.initialize_grid_uniform!(grids)
end

function initialize!(
    model::ModelData,
    ::Val{option_names.TtypeRefinedGrid}
)
    iuse_refinement_delay = model.grids.parameters.refinement.iuse_refinement_delay.value
    if iuse_refinement_delay == 0
        check_t_type_refinement_parameters(model)
        Initialization.initialize_grid_t_type_refined!(model.grids)
    else
        Initialization.initialize_grid_uniform!(model.grids)
    end
end

function initialize!(
    grids::Grids,
    ::Val{option_names.TtypeRefinedGrid}
)
    iuse_refinement_delay = grids.parameters.refinement.iuse_refinement_delay.value
    if iuse_refinement_delay == 0
        check_t_type_refinement_parameters(grids)
        Initialization.initialize_grid_t_type_refined!(grids)
    else
        Initialization.initialize_grid_uniform!(grids)
    end
end

function initialize!(
    grids::Grids3d,
    ::Val{option_names.TtypeRefinedGrid}
)
    iuse_refinement_delay = grids.parameters.refinement.iuse_refinement_delay.value
    if iuse_refinement_delay == 0
        check_t_type_refinement_parameters(grids)
        # TODO: Ttype refinement is not implemented for Grids3d yet
        #Initialization.initialize_grid_t_type_refined!(grids)
        error("Ttype refinement is not implemented for Grids3d yet")
    else
        Initialization.initialize_grid_uniform!(grids)
    end
end

function initialize!(
    model::ModelData,
    ::Val{option_names.Grid51x51UniformRefinementAlongSides}
)
    Initialization.initialize_grid_51x51_uniform_side_refinement!(model.grids)
end

function initialize!(
    model::ModelData,
    ::Val{option_names.Grid51x51SmoothedRefinementAlongSides}
)
    Initialization.initialize_grid_51x51_smoothed_side_refinement!(model.grids)
end

function check_t_type_refinement_parameters(model::ModelData)
    Ttype.GridUpdate.check_t_type_refinement_parameters(model)
end

function check_t_type_refinement_parameters(grids::Grids3d)
    # TODO: Ttype refinement is not implemented for Grids3d yet
    #Ttype.GridUpdate.check_t_type_refinement_parameters(grids)
end

end # module 