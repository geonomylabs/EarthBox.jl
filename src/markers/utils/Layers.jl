module Layers

"""
    calc_depths_for_layers(
        thick_air::Float64,
        thick_upr_crust::Float64,
        thick_lwr_crust::Float64,
        thick_upr_mantle_lithosphere::Float64,
        thick_weak_zone::Float64,
        thick_lwr_mantle_lithosphere::Float64
    ) -> Tuple{Float64, Float64, Float64, Float64, Float64, Float64}

Calculate the base of lithospheric layers.

# Arguments
- `thick_air`: Thickness of the sticky air layer
- `thick_upr_crust`: Thickness of the upper crust
- `thick_lwr_crust`: Thickness of the lower crust
- `thick_upr_mantle_lithosphere`: Thickness of the upper mantle lithosphere
- `thick_weak_zone`: Thickness of the intra-lithosphere weak zone
- `thick_lwr_mantle_lithosphere`: Thickness of the lower mantle lithosphere

# Returns
- `base_air`: Base of the sticky air layer
- `base_upr_crust`: Base of the upper crust
- `base_lwr_crust`: Base of the lower crust
- `base_upr_mantle_lithosphere`: Base of the upper mantle lithosphere
- `base_weak_zone`: Base of the intra-lithosphere weak zone
- `base_lwr_mantle_lithosphere`: Base of the lower mantle lithosphere
"""
function calc_depths_for_layers(
    thick_air::Float64,
    thick_upr_crust::Float64,
    thick_lwr_crust::Float64,
    thick_upr_mantle_lithosphere::Float64,
    thick_weak_zone::Float64,
    thick_lwr_mantle_lithosphere::Float64
)::Tuple{Float64, Float64, Float64, Float64, Float64, Float64}
    base_air = thick_air
    base_upr_crust = base_air + thick_upr_crust
    base_lwr_crust = base_upr_crust + thick_lwr_crust
    base_upr_mantle_lithosphere = base_lwr_crust + thick_upr_mantle_lithosphere
    base_weak_zone = base_upr_mantle_lithosphere + thick_weak_zone
    base_lwr_mantle_lithosphere = base_weak_zone + thick_lwr_mantle_lithosphere
    
    return (
        base_air,
        base_upr_crust,
        base_lwr_crust,
        base_upr_mantle_lithosphere,
        base_weak_zone,
        base_lwr_mantle_lithosphere
    )
end

end # module Layers 