module Regrid

import EarthBox.PrintFuncs: print_info
import EarthBox.ModelDataContainer: ModelData
import EarthBox.PrintFuncs: @timeit_memit
import EarthBox.StaggeredGrid: get_active_option
import EarthBox.StaggeredGrid.Spacing
import EarthBox.StaggeredGrid.Dimensions
import EarthBox.StaggeredGrid.UpdateBasicAndStaggeredGrids
import EarthBox.StaggeredGrid.TypesManager: Uniform
import EarthBox.StaggeredGrid.TypesManager: Ttype

function regrid!(
    model::ModelData,
    boolean_options::Dict{String, Bool}
)::Nothing
    extending_side_boundaries = get(boolean_options, "extending_side_boundaries", false)
    uniform_grid = get(boolean_options, "uniform_grid", false)
    option_active = get_active_option(model)
    if extending_side_boundaries
        advect_grid_boundaries_and_regrid!(model, uniform_grid)
    elseif Symbol(option_active.option_name) == :TtypeRefinedGrid
        @timeit_memit "Finished regridding for TtypeRefinedGrid" begin
            if model.grids.parameters.refinement.iuse_trench.value == 1
                xo_highres_updated = Ttype.TrenchModel.update_highres_region_using_trench!(model)
                print_info("Updated left_edge of highres region using trench model, xo_highres_updated: $xo_highres_updated", level=2)
                Ttype.Regrid.regrid_t_type!(model, xo_highres_updated)
            end

            if model.grids.parameters.refinement.iuse_refinement_delay.value == 1
                Ttype.Regrid.regrid_t_type_with_delay!(model)
            end

            if model.grids.parameters.refinement.iuse_refinement_gap.value == 1
                Ttype.Regrid.regrid_t_type_with_gap!(model)
            end
        end
    end
    return nothing
end

""" Advect grid boundaries, regrid and update bottom velocity bc.

This function only applies to cases with a fixed free slip upper boundary
and moving lower and side boundaries:
    bc_type = 'LithosphericExtensionMovingBoundaries'

No marker recycling is needed for this case.
"""
function advect_grid_boundaries_and_regrid!(model::ModelData, uniform_grid::Bool)
    @timeit_memit "Finished advecting grid boundaries and regridding" begin
        update_model_size!(model)
        Spacing.update_average_spacing_of_basic_grid!(model.grids)
        update_initial_x_location_of_highres_area(model)
        if uniform_grid
            Uniform.Regrid.regrid_uniform_grid_using_extension_velocity!(model)
        else
            Ttype.Regrid.regrid!(model.grids)
        end
        Dimensions.update_grid_dimensions!(model.grids)
        UpdateBasicAndStaggeredGrids.update!(model.grids)
        update_bottom_velocity_bc(model)
    end
end

""" Update model size based on extension velocity.
"""
function update_model_size!(model::ModelData)
    timestep = model.timestep.parameters.main_time_loop.timestep.value
    full_velocity_extension = model.bcs.parameters.velocity.full_velocity_extension.value

    geometry = model.grids.parameters.geometry
    xsize = geometry.xsize.value
    ysize = geometry.ysize.value

    geometry.xsize.value = xsize + full_velocity_extension*timestep
    geometry.ysize.value = ysize - full_velocity_extension/xsize*ysize*timestep
end

""" Advect grid boundaries, regrid and update bottom velocity bc.

This function only applies to cases with a fixed free slip upper boundary
and moving lower and side boundaries:
    bc_type = 'LithosphericExtensionMovingBoundaries'

No marker recycling is needed for this case.
"""
function update_staggered_grid!(model::ModelData, uniform_grid::Bool)
    Spacing.update_average_spacing_of_basic_grid!(model.grids)
    if uniform_grid
        Uniform.Regrid.regrid_uniform_grid_using_extension_velocity!(model)
    else
        Ttype.Regrid.regrid!(model.grids)
    end
    Dimensions.update_grid_dimensions!(model.grids)
    UpdateBasicAndStaggeredGrids.update!(model.grids)
end

""" Update bottom boundary condition for extension model with moving boundaries.

New model width and height requires that the y-component of velocity
is re-calculated. Therefore, the bottom velocity boundary condition needs
to be re-calculated. The boundary condition is free slip with prescribed
inward velocity.

# Updated Arrays
model.bcs.arrays.velocity
    `bbottom.array::Matrix{Float64}` (xnum+1,4):
        Velocity boundary conditions along bottom boundary:
            vx[ynum+1,j] = bbottom[j,1] + vx[ynum,j]*bbottom[j,2]
            vy[ynum+1,j] = bbottom[j,3] + vy[ynum,j]*bbottom[j,4]

model.bcs.arrays.vel_comp
    `bbottomy.array::Matrix{Float64}` (xnum+1,2):
        Velocity boundary conditions along bottom boundary:
            vy[ynum+1,j] = bbottomy[j,1] + vy[ynum,j]*bbottomy[j,2]
"""
function update_bottom_velocity_bc(model::ModelData)::Nothing
    xsize = model.grids.parameters.geometry.xsize.value
    ysize = model.grids.parameters.geometry.ysize.value
    xnum = model.grids.parameters.geometry.xnum.value
    full_velocity_extension = model.bcs.parameters.velocity.full_velocity_extension.value

    bbottom = model.bcs.arrays.velocity.bbottom.array
    bbottomy = model.bcs.arrays.vel_comp.bbottomy.array

    for j in 1:xnum+1
        bbottom[j,3] = -full_velocity_extension/xsize*ysize
        bbottom[j,4] = 0.0
        bbottomy[j,1] = -full_velocity_extension/xsize*ysize
        bbottomy[j,2] = 0.0
    end

    return nothing
end

end # module

