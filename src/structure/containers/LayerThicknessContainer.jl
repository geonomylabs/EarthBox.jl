module LayerThicknessContainer

mutable struct LayerThickness
    # Thickness of sticky layer in meters
    thick_sticky::Union{Vector{Float64}, Nothing}
    # Thickness of crust in meters
    thick_crust::Union{Vector{Float64}, Nothing}
    # Thickness of mantle lithosphere in meters
    thick_mantle_lith::Union{Vector{Float64}, Nothing}
    # Thickness of asthenosphere in meters
    thick_asthenosphere::Union{Vector{Float64}, Nothing}
    # Thickness of sediments in meters
    thick_sediments::Union{Vector{Float64}, Nothing}
end

# Constructor with default values
LayerThickness(; thick_sticky=nothing, thick_crust=nothing,
    thick_mantle_lith=nothing, thick_asthenosphere=nothing,
    thick_sediments=nothing) = 
    LayerThickness(thick_sticky, thick_crust, thick_mantle_lith,
        thick_asthenosphere, thick_sediments)

end # module 