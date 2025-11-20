""" Print information about the lava flow model.
"""
module PrintLavaFlowInfo

import EarthBox.PrintFuncs: print_flow_info
export print_lava_model_info, print_flow_info

"""
    print_lava_model_info(
        idrainage_basin::Int,
        xmid_molten_domain::Float64,
        width_eruption_domain::Float64,
        eruption_location_x_min::Float64,
        eruption_location_x_forecast::Float64,
        eruption_location_y_forecast::Float64,
        total_extrusion_volume::Float64,
        characteristic_volume_per_flow::Float64,
        number_of_flows_per_model_time_step::Int
    )

Print information about the lava flow model.
"""
function print_lava_model_info(
    idrainage_basin::Int,
    xmid_molten_domain::Float64,
    width_eruption_domain::Float64,
    eruption_location_x_min::Float64,
    eruption_location_x_forecast::Float64,
    eruption_location_y_forecast::Float64,
    total_extrusion_volume::Float64,
    characteristic_volume_per_flow::Float64,
    number_of_flows_per_model_time_step::Int
)
    print_flow_info("", level=1)
    print_flow_info("Working on flows for drainage basin $idrainage_basin", level=2)
    print_flow_info("xmid_molten_domain: $xmid_molten_domain", level=3)  
    print_flow_info("width_eruption_domain: $width_eruption_domain", level=3)
    print_flow_info("eruption_location_x_min: $eruption_location_x_min", level=3)
    print_flow_info("eruption_location_x_forecast: $eruption_location_x_forecast", level=3)
    print_flow_info("eruption_location_y_forecast: $eruption_location_y_forecast", level=3)
    print_flow_info("total extrusion volume: $total_extrusion_volume", level=3)
    print_flow_info("characteristic volume per flow: $characteristic_volume_per_flow", level=3)
    print_flow_info("number of flows per model time step: $number_of_flows_per_model_time_step", level=3)
    print_flow_info("", level=1)
end

"""
    print_flow_info(
        mm::Int,
        eruption_location_x::Float64,
        eruption_location_y::Float64,
        y_sealevel::Float64,
        eruption_stype::String,
        residual_lava_thickness::Float64,
        flow_volume::Float64,
        flow_thickness_max::Float64
    )

Print information about the lava flow.
"""
function print_flow_info(
    mm::Int,
    eruption_location_x::Float64,
    eruption_location_y::Float64,
    y_sealevel::Float64,
    eruption_stype::String,
    residual_lava_thickness::Float64,
    flow_volume::Float64,
    flow_thickness_max::Float64
)
    print_flow_info("", level=2)
    print_flow_info("Flow number: $mm", level=3)
    print_flow_info("Eruption x-coordinate: $eruption_location_x", level=3)
    print_flow_info("Eruption y-coordinate: $eruption_location_y", level=3)
    print_flow_info("Eruption y_sealevel: $y_sealevel", level=3)
    print_flow_info("Eruption type: $eruption_stype", level=3)
    print_flow_info("Flow volume: $flow_volume", level=3)
    print_flow_info("Residual lava thickness: $residual_lava_thickness", level=3)
    print_flow_info("Maximum flow calculated thickness: $flow_thickness_max", level=3)
end

end # module 