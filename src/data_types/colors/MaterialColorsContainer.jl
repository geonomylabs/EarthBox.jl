module MaterialColorsContainer

"""
    MaterialColors

Struct with material color information.

# Fields
- `color_map::`[`Union{Any, Nothing}`]:
    - Color map for material colors
- `n_bin::`[`Union{Int64, Nothing}`]:
    - Number of material colors
- `labels::`[`Union{Vector{String}, Nothing}`]:
    - Labels for material colors with format: `material type : material domain : (material name)`
- `ticks::`[`Union{Vector{Float64}, Nothing}`]:
    - Ticks for the colorbar
- `colorrange::`[`Union{Tuple{Float64, Float64}, Nothing}`]: 
    - Color range for the colorbar specifically for CairoMakie where for material
       ID's ranging from `1`` to `n_bin` the colorbar range should be from `1 - 0.5`
       to `n_bin + 0.5` to maintain proper spacing between the colorbar ticks.
"""
mutable struct MaterialColors
    color_map::Union{Any, Nothing}
    n_bin::Union{Int64, Nothing}
    labels::Union{Vector{String}, Nothing}
    ticks::Union{Vector{Float64}, Nothing}
    colorrange::Union{Tuple{Float64, Float64}, Nothing}
end

function MaterialColors()
    return MaterialColors(nothing, nothing, nothing, nothing, nothing)
end

function get_colorbar_ticks_for_material_colors(n_bin::Int)::Vector{Float64}
    ticks = Float64[]
    for i in 1:n_bin
        push!(ticks, Float64(i))
    end
    return ticks
end


end