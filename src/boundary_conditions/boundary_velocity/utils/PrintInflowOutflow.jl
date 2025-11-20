module PrintInflowOutflow

import EarthBox.ConversionFuncs: meters_per_seconds_to_cm_per_yr

function print_inflow_velocity(
    velocity_y_upper::Float64,
    velocity_y_lower::Float64
)
    println(
        ">> velocity_y_upper (cm/yr): ",
        meters_per_seconds_to_cm_per_yr(velocity_y_upper)
    )
    println(
        ">> velocity_y_lower (cm/yr): ",
        meters_per_seconds_to_cm_per_yr(velocity_y_lower)
    )
end

function print_average_side_velocity(
    velocity_left_avg::Float64,
    velocity_right_avg::Float64
)
    println(
        ">> velocity_left_avg (cm/yr): ",
        meters_per_seconds_to_cm_per_yr(velocity_left_avg)
    )
    println(
        ">> velocity_right_avg (cm/yr): ",
        meters_per_seconds_to_cm_per_yr(velocity_right_avg)
    )
end

end # module PrintInflowOutflow 