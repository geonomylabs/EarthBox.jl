module Regrid

import EarthBox.ModelDataContainer.Grids2dContainer: Grids
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ConversionFuncs: seconds_to_years
import EarthBox.ParameterGroupTools: set_group_parameters!
import ...Uniform.Initialization: initialize! as initialize_grid_uniform!
import ..Ttype.Initialization: get_left_edge_x_index_of_highres_zone
import ..Ttype.Initialization: calculate_final_xgrid_index_of_highres
import ..Ttype.Initialization: define_xgrid_to_the_right_of_high_res_region!
import ..Ttype.Initialization: define_xgrid_to_the_left_of_high_res_region!
import ..Ttype.Initialization: calculate_gridy_coordinates!
import ..Ttype.GridUpdate: update_t_type_grid!

function regrid!(grids::Grids)
    regrid_gridx!(grids)
    calculate_gridy_coordinates!(grids)
end

function regrid_gridx!(grids::Grids)
    # Initial x-grid location of high-resolution area
    xo_highres = grids.parameters.refinement.xo_highres.value
    # Final x-grid location of high-resolution area
    xf_highres = grids.parameters.refinement.xf_highres.value
    # X-direction spacing in high-resolution area
    dx_highres = grids.parameters.refinement.dx_highres.value
    # Number of nodes of basic grid in x-direction
    xnum = grids.parameters.geometry.xnum.value
    # Size of model in x-direction
    xsize = grids.parameters.geometry.xsize.value
    # X-coordinates of basic grid in x-direction
    gridx_b = grids.arrays.basic.gridx_b.array

    # Initial x-grid index of high-resolution area
    ixo_highres = get_left_edge_x_index_of_highres_zone(
        xnum, xo_highres, xf_highres, dx_highres)
    # Final x-grid index of high-resolution area
    ixf_highres = calculate_final_xgrid_index_of_highres(
        ixo_highres, xo_highres, xf_highres, dx_highres)

    # Define position of first nodal point
    # The following update for 1-based indexing has not been tested
    ixm = Int64(floor((xnum-1)/2)) + 1 # +1 to account for 1-based indexing
    gridx_b[1] = gridx_b[ixm] - xsize/2.0

    # Distance to be covered by right-side non-uniform grid
    dist = xsize/2.0 - (gridx_b[ixf_highres] - gridx_b[ixm])
    define_xgrid_to_the_right_of_high_res_region!(
        xnum, ixf_highres, dist, dx_highres, gridx_b)

    # Distance to be covered by left-side non-uniform grid
    dist = gridx_b[ixo_highres] - gridx_b[1]
    define_xgrid_to_the_left_of_high_res_region!(
        ixo_highres, dist, dx_highres, gridx_b)
end

""" Update initial x-location of high-resolution area.

This function applies to cases with a moving left-side boundary.
"""
function update_initial_x_location_of_highres_area(model::ModelData)
    timestep = model.timestep.parameters.main_time_loop.timestep.value
    full_velocity_extension = model.bcs.parameters.velocity.full_velocity_extension.value
    gridx_b = model.grids.arrays.basic.gridx_b.array

    xo_initial = gridx_b[1]
    xo_final = xo_initial - full_velocity_extension/2.0*timestep

    model.grids.parameters.refinement.xo_highres.value = xo_final
end

function regrid_t_type!(model::ModelData, xo_highres_updated::Float64)
    xo_highres = model.grids.parameters.refinement.xo_highres.value
    xf_highres = model.grids.parameters.refinement.xf_highres.value
    width_highres = xf_highres - xo_highres
    xf_highres_updated = xo_highres_updated + width_highres

    set_group_parameters!(
        model.grids.parameters.refinement, 
        Dict(
            "xo_highres" => xo_highres_updated,
            "xf_highres" => xf_highres_updated
            )
        )
    update_t_type_grid!(model)
end

function regrid_t_type_with_delay!(model::ModelData)
    timesum = model.timestep.parameters.main_time_loop.timesum.value
    timesum_myr = seconds_to_years(timesum) / 1e6
    refinement_time_myr = model.grids.parameters.refinement.refinement_time.value
    refinement_flag = model.grids.parameters.refinement.refinement_flag.value

    if timesum_myr >= refinement_time_myr && refinement_flag == 0
        println(">> Refining grid using t-type refinement.")
        update_t_type_grid!(model)
        update_refinement_flag(model, 1)
    end
end

function update_refinement_flag(model::ModelData, flag::Int)
    model.grids.parameters.refinement.refinement_flag.value = flag
end

function regrid_t_type_with_gap!(model::ModelData)
    timesum = model.timestep.parameters.main_time_loop.timesum.value
    timesum_myr = seconds_to_years(timesum) / 1e6
    gap_start = model.grids.parameters.refinement.refinement_gap_start_time.value
    gap_end = model.grids.parameters.refinement.refinement_gap_end_time.value

    if gap_start <= timesum_myr <= gap_end
        println(">> Regridding using uniform grid within refinement gap.")
        initialize_grid_uniform!(model.grids)
    else
        println(">> Regridding using t-type refinement.")
        update_t_type_grid!(model)
    end
end

end # module TtypeRegrid

