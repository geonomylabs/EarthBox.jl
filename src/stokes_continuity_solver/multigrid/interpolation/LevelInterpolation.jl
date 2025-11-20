module LevelInterpolation

import ..MultigridDataManager: MultigridData3d, MultigridData2d
import ..ViscosityRestriction: viscosity_restriction3d!, viscosity_restriction2d!

function interpolate_viscosity_to_coarser_levels!(
    multigrid_data::MultigridData3d
)::Nothing
    level_vector = multigrid_data.level_vector
    levelnum = length(level_vector)
    # Interpolation (restriction) of viscosity to coarser levels
    for n = 1:levelnum-1
        viscosity_restriction3d!(n, level_vector)
    end
    return nothing
end

function interpolate_viscosity_to_coarser_levels!(
    multigrid_data::MultigridData2d
)::Nothing
    level_vector = multigrid_data.level_vector
    levelnum = length(level_vector)
    # Interpolation (restriction) of viscosity to coarser levels
    for n = 1:levelnum-1
        viscosity_restriction2d!(n, level_vector)
    end
    return nothing
end

end # module