module ConservationCheck

import EarthBox.ModelDataContainer: ModelData
import ..AverageVelocity: calculate_average_sub_plate_depth_dependent_velocity
import ..OutflowGeometry: get_outflow_geometry, get_outflow_size
import ..InflowGeometry: get_inflow_geometry
import ..StickyThickness: get_sticky_thickness

"""
    check_conservation_constant(
        xsize::Float64, ysize::Float64, velocity_y_upper::Float64, 
        velocity_y_lower::Float64, velocity_left::Float64, 
        velocity_right::Float64, sticky_thickness_left::Float64, 
        sticky_thickness_right::Float64
    )::Nothing

Check mass conservation at top and bottom boundaries of extension model.

This function applies only to cases with uniform extension along side boundaries.
"""
function check_conservation_constant(
    xsize::Float64,
    ysize::Float64,
    velocity_y_upper::Float64,
    velocity_y_lower::Float64,
    velocity_left::Float64,
    velocity_right::Float64,
    sticky_thickness_left::Float64,
    sticky_thickness_right::Float64
)
    inflow = calc_inflow(xsize, velocity_y_upper, velocity_y_lower)

    outflow = calc_outflow(ysize, ysize, velocity_left, velocity_right)

    sticky_outflow = calc_outflow(
        sticky_thickness_left, sticky_thickness_right,
        velocity_left, velocity_right
    )

    rock_outflow = calculate_rock_outflow_constant(
        ysize, velocity_left, velocity_right,
        sticky_thickness_left, sticky_thickness_right
    )

    outflow_sticky_rock = sticky_outflow + rock_outflow
    net_flow = inflow - outflow

    print_conservation_warnings(net_flow, outflow, outflow_sticky_rock)
end

"""
    calculate_rock_outflow_constant(
        ysize::Float64, velocity_left::Float64, velocity_right::Float64, 
        sticky_thickness_left::Float64, sticky_thickness_right::Float64
    )::Float64

Calculate rock outflow at top and bottom boundaries of extension model.

This function applies only to cases with uniform extension along side boundaries.
"""
function calculate_rock_outflow_constant(
    ysize::Float64,
    velocity_left::Float64,
    velocity_right::Float64,
    sticky_thickness_left::Float64,
    sticky_thickness_right::Float64
)::Float64
    rock_thickness_left = ysize - sticky_thickness_left
    rock_thickness_right = ysize - sticky_thickness_right

    rock_outflow = calc_outflow(
        rock_thickness_left, rock_thickness_right,
        velocity_left, velocity_right
    )

    return rock_outflow
end

"""
    check_conservation_depth_dependent(
        xsize::Float64, ysize::Float64, y_linear_left::Float64, y_linear_right::Float64, 
        velocity_avg_left::Float64, velocity_avg_right::Float64, velocity_y_upper::Float64, 
        velocity_y_lower::Float64, velocity_left::Float64, velocity_right::Float64, 
        sticky_thickness_left::Float64, sticky_thickness_right::Float64
    )::Nothing

Check mass conservation at top and bottom boundaries of extension model.

This function applies only to cases with depth-dependent extension along side boundaries.
"""
function check_conservation_depth_dependent(
    xsize::Float64,
    ysize::Float64,
    y_linear_left::Float64,
    y_linear_right::Float64,
    velocity_avg_left::Float64,
    velocity_avg_right::Float64,
    velocity_y_upper::Float64,
    velocity_y_lower::Float64,
    velocity_left::Float64,
    velocity_right::Float64,
    sticky_thickness_left::Float64,
    sticky_thickness_right::Float64
)
    inflow = calc_inflow(xsize, velocity_y_upper, velocity_y_lower)

    outflow = calc_outflow(ysize, ysize, velocity_avg_left, velocity_avg_right)

    sticky_outflow = calc_outflow(
        sticky_thickness_left, sticky_thickness_right,
        velocity_left, velocity_right
    )

    rock_outflow = calculate_rock_outflow_depth_dependent(
        ysize, y_linear_left, y_linear_right, velocity_left, velocity_right,
        sticky_thickness_left, sticky_thickness_right
    )

    outflow_sticky_rock = sticky_outflow + rock_outflow
    net_flow = inflow - outflow

    print_conservation_warnings(net_flow, outflow, outflow_sticky_rock)
end

"""
    calculate_rock_outflow_depth_dependent(
        ysize::Float64, y_linear_left::Float64, y_linear_right::Float64, 
        velocity_left::Float64, velocity_right::Float64, 
        sticky_thickness_left::Float64, sticky_thickness_right::Float64
    )::Float64

Calculate rock outflow at top and bottom boundaries of extension model.

This function applies only to cases with depth-dependent extension along side boundaries.
"""
function calculate_rock_outflow_depth_dependent(
    ysize::Float64,
    y_linear_left::Float64,
    y_linear_right::Float64,
    velocity_left::Float64,
    velocity_right::Float64,
    sticky_thickness_left::Float64,
    sticky_thickness_right::Float64
)::Float64
    plate_thickness_left = y_linear_left - sticky_thickness_left
    plate_thickness_right = y_linear_right - sticky_thickness_right

    plate_rock_outflow = calc_outflow(
        plate_thickness_left, plate_thickness_right,
        velocity_left, velocity_right
    )

    (
        velocity_avg_below_plate_left, 
        velocity_avg_below_plate_right
    ) = calc_subplate_velocity(velocity_left, velocity_right)

    subplate_rock_outflow = calc_outflow(
        ysize-sticky_thickness_left-plate_thickness_left,
        ysize-sticky_thickness_right-plate_thickness_right,
        velocity_avg_below_plate_left,
        velocity_avg_below_plate_right
    )

    return plate_rock_outflow + subplate_rock_outflow
end

"""
    calc_subplate_velocity(
        velocity_left::Float64, velocity_right::Float64)::Tuple{Float64, Float64}

Calculate average subplate velocity.
"""
function calc_subplate_velocity(
    velocity_left::Float64,
    velocity_right::Float64
)::Tuple{Float64, Float64}
    (
        velocity_avg_below_plate_left
    ) = calculate_average_sub_plate_depth_dependent_velocity(velocity_left)
    (
        velocity_avg_below_plate_right
    ) = calculate_average_sub_plate_depth_dependent_velocity(velocity_right)
    return (velocity_avg_below_plate_left, velocity_avg_below_plate_right)
end

"""
    calc_inflow(
        xsize::Float64, velocity_y_upper::Float64, velocity_y_lower::Float64
    )::Float64

Calculate inflow along top and bottom boundaries.

inflow = velocity_y_upper*xsize - velocity_y_lower*xsize
"""
function calc_inflow(
    xsize::Float64,
    velocity_y_upper::Float64,
    velocity_y_lower::Float64
)::Float64
    return velocity_y_upper*xsize - velocity_y_lower*xsize
end

"""
    calc_outflow(
        ysize_left::Float64, ysize_right::Float64, velocity_left::Float64, 
        velocity_right::Float64
    )::Float64

Calculate outflow along left and right boundaries.

outflow = velocity_right*ysize_right - velocity_left*ysize_left
"""
function calc_outflow(
    ysize_left::Float64,
    ysize_right::Float64,
    velocity_left::Float64,
    velocity_right::Float64
)::Float64
    return velocity_right*ysize_right - velocity_left*ysize_left
end

"""
    print_conservation_warnings(
        net_flow::Float64, outflow::Float64, outflow_sticky_rock::Float64
    )::Nothing

Print conservation warnings if applicable.

# Arguments
- `net_flow::Float64`: 
    - Net flow of model velocity bcs: inflow - outflow
- `outflow::Float64`: 
    - Outflow along model boundaries using constant or average side boundary velocity
- `outflow_sticky_rock::Float64`: 
    - Outflow along sticky and rock boundaries calculated taking into account 
      each boundary segment separately
"""
function print_conservation_warnings(
    net_flow::Float64,
    outflow::Float64,
    outflow_sticky_rock::Float64
)
    if abs(net_flow) > 1e-6
        println(
            ">> !!! WARNING !!! Net flow of velocity boundary conditions is " *
            "not zero: ", net_flow
        )
    end

    if abs(outflow-outflow_sticky_rock) > 1e-6
        println(
            ">> !!! WARNING !!! Outflow calculated using average side " *
            "velocity boundary conditions does not match outflow calculated " *
            "using separate segments along sticky and rock boundaries: ",
            outflow, " ", outflow_sticky_rock
        )
    end
end

"""
    check_conservation_inflow_and_outflow_along_sides(
        model::ModelData, velocity_y_upper::Float64, 
        velocity_avg_outflow_left::Float64, velocity_avg_outflow_right::Float64, 
        velocity_left::Float64, velocity_right::Float64
    )::Nothing

Check mass conservation at top and bottom boundaries of extension model.
"""
function check_conservation_inflow_and_outflow_along_sides(
    model::ModelData,
    velocity_y_upper::Float64,
    velocity_avg_outflow_left::Float64,
    velocity_avg_outflow_right::Float64,
    velocity_left::Float64,
    velocity_right::Float64
)::Nothing
    xsize = model.grids.parameters.geometry.xsize.value
    plate_thickness, smoothing_thickness = get_outflow_geometry(model)
    outflow_size_left, outflow_size_right = get_outflow_size(model)
    ysize_inflow_left, ysize_inflow_right = get_inflow_geometry(model)
    sticky_thickness_left, sticky_thickness_right = get_sticky_thickness(model)
    velocity_inflow_left = model.bcs.parameters.velocity.velocity_inflow_left.value
    velocity_inflow_right = model.bcs.parameters.velocity.velocity_inflow_right.value

    inflow = calc_inflow_top_and_sides(
        xsize, velocity_y_upper, smoothing_thickness,
        ysize_inflow_left, ysize_inflow_right,
        velocity_inflow_left, velocity_inflow_right
        )

    outflow = calc_outflow_sides(
        velocity_avg_outflow_left,
        velocity_avg_outflow_right,
        outflow_size_left,
        outflow_size_right
        )

    sticky_outflow = calc_outflow(
        sticky_thickness_left, sticky_thickness_right,
        velocity_left, velocity_right
    )

    rock_outflow = calculate_rock_outflow_sides(
        plate_thickness, smoothing_thickness, velocity_left, velocity_right)

    outflow_sticky_rock = sticky_outflow + rock_outflow
    net_flow = inflow - outflow

    print_conservation_warnings(net_flow, outflow, outflow_sticky_rock)
    return nothing
end

"""
    calc_inflow_top_and_sides(
        xsize::Float64, velocity_y_upper::Float64, smoothing_thickness::Float64, 
        ysize_inflow_left::Float64, ysize_inflow_right::Float64, 
        velocity_inflow_left::Float64, velocity_inflow_right::Float64
    )::Float64

Calculate inflow along top and sides.
"""
function calc_inflow_top_and_sides(
    xsize::Float64,
    velocity_y_upper::Float64,
    smoothing_thickness::Float64,
    ysize_inflow_left::Float64,
    ysize_inflow_right::Float64,
    velocity_inflow_left::Float64,
    velocity_inflow_right::Float64
)::Float64
    ysize_inflow_left_constant = ysize_inflow_left - smoothing_thickness
    ysize_inflow_right_constant = ysize_inflow_right - smoothing_thickness

    inflow = (
        velocity_y_upper*xsize
        + 0.5*velocity_inflow_left*smoothing_thickness
        + velocity_inflow_left*ysize_inflow_left_constant
        - 0.5*velocity_inflow_right*smoothing_thickness
        - velocity_inflow_right*ysize_inflow_right_constant
        )
    return inflow
end

"""
    calc_outflow_sides(
        velocity_avg_outflow_left::Float64, velocity_avg_outflow_right::Float64, 
        outflow_size_left::Float64, outflow_size_right::Float64
    )::Float64

Calculate outflow along sides.
"""
function calc_outflow_sides(
    velocity_avg_outflow_left::Float64,
    velocity_avg_outflow_right::Float64,
    outflow_size_left::Float64,
    outflow_size_right::Float64
)::Float64
    return (
        velocity_avg_outflow_right*outflow_size_right 
        - velocity_avg_outflow_left*outflow_size_left
        )
end

"""
    calculate_rock_outflow_sides(
        plate_thickness::Float64, smoothing_thickness::Float64, 
        velocity_left::Float64, velocity_right::Float64
    )::Float64

Calculate rock outflow along sides.
"""
function calculate_rock_outflow_sides(
    plate_thickness::Float64,
    smoothing_thickness::Float64,
    velocity_left::Float64,
    velocity_right::Float64
)::Float64
    plate_rock_outflow = calc_outflow(
        plate_thickness, plate_thickness,
        velocity_left, velocity_right
        )

    subplate_rock_outflow = calc_outflow(
        smoothing_thickness, smoothing_thickness,
        velocity_left*0.5,
        velocity_right*0.5
        )

    return plate_rock_outflow + subplate_rock_outflow
end

end # module ConservationCheck 