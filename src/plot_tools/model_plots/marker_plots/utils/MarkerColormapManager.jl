module MarkerColormapManager

import CairoMakie
import EarthBox.Markers.MarkerMaterials.MaterialsContainer: Materials
import EarthBox.Markers.MarkerMaterials.MaterialsContainer.MaterialsStateContainer: MaterialsState
import EarthBox.Markers.MarkerMaterials.MaterialsContainer.MaterialsStateContainer: isloaded
import ..Ticks: get_colorbar_ticks_for_composition_plot


mutable struct MarkerColorMap
    cm::Union{Any, Nothing}
    n_bin::Int64
end

function MarkerColorMap(
    materials::Materials,
    materials_state::MaterialsState
)::MarkerColorMap
    cm = nothing
    n_bin = 0
    if isloaded(materials_state)
        cm, n_bin = make_color_map_comp(materials)
    end
    return MarkerColorMap(cm, n_bin)
end

function make_color_map_comp(materials::Materials)::Tuple{Any, Int64}
    materials_dict = materials.materials
    material_ids = collect(keys(materials_dict))
    ncolors = maximum(material_ids)
    colors = Vector{CairoMakie.RGB{Float64}}()
    
    for matid in 1:ncolors
        if matid in material_ids
            rgb = materials_dict[matid].rgb
            red_fraction = rgb.red_fraction.value
            green_fraction = rgb.green_fraction.value
            blue_fraction = rgb.blue_fraction.value
        else
            red_fraction = 0.0
            green_fraction = 0.0
            blue_fraction = 0.0
        end
        check_rbs_color_values(matid, red_fraction, green_fraction, blue_fraction)
        color = CairoMakie.RGB(red_fraction, green_fraction, blue_fraction)
        push!(colors, color)
    end
    
    n_bin = ncolors
    if length(colors) < 2
        push!(colors, CairoMakie.RGB(1.0, 1.0, 1.0))
        n_bin += 1
    end
    
    cm = CairoMakie.cgrad(colors; categorical = true)
    return cm, n_bin
end

function check_rbs_color_values(
    matid::Int64,
    red_fraction::Float64, 
    green_fraction::Float64, 
    blue_fraction::Float64
)::Nothing
    # Ensure all values are Float64
    if !isa(red_fraction, Float64) || !isa(green_fraction, Float64) || !isa(blue_fraction, Float64)
        error("RGB fractions must be Float64. Check material input for material $matid.")
    end
    # Check for NaN and valid range 0 to 1
    if isnan(red_fraction) || isnan(green_fraction) || isnan(blue_fraction)
        error("RGB fractions must not be NaN. Check material input for material $matid.")
    end
    if red_fraction < 0.0 || red_fraction > 1.0 ||
       green_fraction < 0.0 || green_fraction > 1.0 ||
       blue_fraction < 0.0 || blue_fraction > 1.0
        error("RGB fractions must be floats between 0 and 1. Check material input for material $matid.")
    end
    return nothing
end

end # module
