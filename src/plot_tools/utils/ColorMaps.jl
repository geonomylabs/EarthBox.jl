module ColorMaps

import CairoMakie

function make_transparent_single_color_colormap(
    color::Tuple{Float64, Float64, Float64};
    max_alpha::Float64=1.0
)::Any
    colors = Vector{CairoMakie.RGBA{Float64}}()
    
    for i in 1:100
        alpha = (i - 1) * max_alpha / 99.0  # 0-based to 1-based conversion
        rgba = CairoMakie.RGBA(color[1], color[2], color[3], alpha)
        push!(colors, rgba)
    end
    
    return CairoMakie.cgrad(colors)
end

function make_alternating_colormap(
    value_min::Float64,
    value_max::Float64,
    delta_value::Float64,
    color1::Tuple{Float64, Float64, Float64},
    color2::Tuple{Float64, Float64, Float64}
)::Any
    num_alternations = floor(Int, (value_max - value_min) / (2.0 * delta_value))
    colors = Vector{CairoMakie.RGB{Float64}}()
    
    for i in 1:num_alternations
        push!(colors, CairoMakie.RGB(color1[1], color1[2], color1[3]))
        push!(colors, CairoMakie.RGB(color2[1], color2[2], color2[3]))
    end
    
    return CairoMakie.cgrad(colors; categorical=true)
end

function make_discontinuous_colormap_from_input_cmap(
    value_min::Float64,
    value_max::Float64,
    delta_value::Float64,
    cmap::Any
)::Any
    # Convert string or symbol to actual colormap if needed
    colormap = if isa(cmap, String)
        CairoMakie.cgrad(Symbol(cmap))
    elseif isa(cmap, Symbol)
        CairoMakie.cgrad(cmap)
    else
        cmap
    end
    # Calculate the number of discrete color steps
    num_steps = floor(Int, (value_max - value_min) / delta_value) + 1
    # Sample colors from the input colormap at regular intervals
    colors = Vector{CairoMakie.RGB{Float64}}()
    for i in 0:(num_steps-1)
        # Sample position in the colormap (0.0 to 1.0)
        t = num_steps > 1 ? i / (num_steps - 1) : 0.5
        # Get the color at this position from the input colormap
        sampled_color = CairoMakie.get(colormap, t)
        push!(colors, CairoMakie.RGB(sampled_color.r, sampled_color.g, sampled_color.b))
    end
    return CairoMakie.cgrad(colors; categorical=true)
end

"""
    struct CairoMakieColorMapNames

Registry of CairoMakie color map names (This needs to be updated for CairoMakie).

Descriptions are from:
https://www.adamsmith.haus/python/docs/matplotlib.pyplot.colormaps

Colormaps derived from those of the same name provided with Matlab:

- autumn: sequential linearly-increasing shades of red-orange-yellow
- bone: sequential increasing black-white color map with a tinge of blue, to 
        emulate X-ray film
- cool: linearly-decreasing shades of cyan-magenta
- copper: sequential increasing shades of black-copper
- flag: repetitive red-white-blue-black pattern (not cyclic at endpoints)
- gray: sequential linearly-increasing black-to-white grayscale
- hot: sequential black-red-yellow-white, to emulate blackbody radiation from 
       an object at increasing temperatures
- hsv: cyclic red-yellow-green-cyan-blue-magenta-red, formed by changing the 
       hue component in the HSV color space
- inferno: perceptually uniform shades of black-red-yellow
- jet: a spectral map with dark endpoints, blue-cyan-yellow-red; based on a 
       fluid-jet simulation by NCSA
- magma: perceptually uniform shades of black-red-white
- pink: sequential increasing pastel black-pink-white, meant for sepia tone 
        colorization of photographs
- plasma: perceptually uniform shades of blue-red-yellow
- prism: repetitive red-yellow-green-blue-purple-...-green pattern (not cyclic 
         at endpoints)
- spring: linearly-increasing shades of magenta-yellow
- summer: sequential linearly-increasing shades of green-yellow
- viridis: perceptually uniform shades of blue-green-yellow
- winter: linearly-increasing shades of blue-green

Palettes from the Yorick scientific visualization package, an evolution of
the GIST package, both by David H. Munro:

- gist_earth: mapmaker's colors from dark blue deep ocean to green lowlands to 
              brown highlands to white mountains
- gist_heat: sequential increasing black-red-orange-white, to emulate 
             blackbody radiation from an iron bar as it grows hotter
- gist_ncar: pseudo-spectral black-blue-green-yellow-red-purple-white colormap 
             from National Center for Atmospheric Research
- gist_rainbow: runs through the colors in spectral order from red to violet at 
                full saturation (like hsv but not cyclic)
- gist_stern: "Stern special" color table from Interactive Data Language 
              software

Colormaps based on the ColorBrewer color specifications and designed and
developed by Cynthia Brewer:

- rBG: brown, white, blue-green
- PiYG: pink, white, yellow-green
- PRGn: purple, white, green
- PuOr: orange, white, purple
- RdBu: red, white, blue
- RdGy: red, white, gray
- RdYlBu: red, yellow, blue
- RdYlGn: red, yellow, green
- Spectral: red, orange, yellow, green, blue

ColorBrewer Sequential (luminance decreases monotonically):

- Blues: white to dark blue
- BuGn: white, light blue, dark green
- BuPu: white, light blue, dark purple
- GnBu: white, light green, dark blue
- Greens: white to dark green
- Greys: white to black (not linear)
- Oranges: white, orange, dark brown
- OrRd: white, orange, dark red
- PuBu: white, light purple, dark blue
- PuBuGn: white, light purple, dark green
- PuRd: white, light purple, dark red
- Purples: white to dark purple
- RdPu: white, pink, dark purple
- Reds: white to dark red
- YlGn: light yellow, dark green
- YlGnBu: light yellow, light green, dark blue
- YlOrBr: light yellow, orange, dark brown
- YlOrRd: light yellow, orange, dark red

Other miscellaneous schemes:

- afmhot: sequential black-orange-yellow-white blackbody spectrum, commonly 
          used in atomic force microscopy
- brg: blue-red-green
- bwr: diverging blue-white-red
- coolwarm: diverging blue-gray-red, meant to avoid issues with 3D shading, 
            color blindness, and ordering of colors
- CMRmap: "Default colormaps on color images often reproduce to confusing 
          grayscale images. The proposed colormap maintains an aesthetically 
          pleasing color image that automatically reproduces to a monotonic 
          grayscale with discrete, quantifiable saturation levels."
- cubehelix: Unlike most other color schemes cubehelix was designed by D.A. 
             Green to be monotonically increasing in terms of perceived 
             brightness. Also, when printed on a black and white postscript 
             printer, the scheme results in a greyscale with monotonically 
             increasing brightness. This color scheme is named cubehelix 
             because the r,g,b values produced can be visualised as a squashed 
             helix around the diagonal in the r,g,b color cube.
- gnuplot: gnuplot's traditional pm3d scheme (black-blue-red-yellow)
- gnuplot2: sequential color printable as gray (black-blue-violet-yellow-white)
- ocean: green-blue-white
- rainbow: spectral purple-blue-green-yellow-orange-red colormap with 
           diverging luminance
- seismic: diverging blue-white-red
- nipy_spectral: black-purple-blue-green-yellow-red-white spectrum, originally 
                 from the Neuroimaging in Python project
- terrain: mapmaker's colors, blue-green-yellow-brown-white, originally from 
           IGOR Pro
"""
struct CairoMakieColorMapNames
    # Matplotlib-derived colormaps
    autumn::String
    bone::String
    cool::String
    copper::String
    flag::String
    gray::String
    hot::String
    hsv::String
    inferno::String
    jet::String
    magma::String
    pink::String
    plasma::String
    prism::String
    spring::String
    summer::String
    viridis::String
    winter::String
    
    # GIST colormaps
    gist_earth::String
    gist_heat::String
    gist_ncar::String
    gist_rainbow::String
    gist_stern::String
    
    # ColorBrewer diverging
    rBG::String
    PiYG::String
    PRGn::String
    PuOr::String
    RdBu::String
    RdGy::String
    RdYlBu::String
    RdYlGn::String
    Spectral::String
    
    # ColorBrewer sequential
    Blues::String
    BuGn::String
    BuPu::String
    GnBu::String
    Greens::String
    Greys::String
    Oranges::String
    OrRd::String
    PuBu::String
    PuBuGn::String
    PuRd::String
    Purples::String
    RdPu::String
    Reds::String
    YlGn::String
    YlGnBu::String
    YlOrBr::String
    YlOrRd::String
    
    # Miscellaneous
    afmhot::String
    brg::String
    bwr::String
    coolwarm::String
    CMRmap::String
    cubehelix::String
    gnuplot::String
    gnuplot2::String
    ocean::String
    rainbow::String
    seismic::String
    nipy_spectral::String
    terrain::String
end

# Default instance with all color map names
const COLORMAP_NAMES = CairoMakieColorMapNames(
    # Matplotlib-derived colormaps
    "autumn", "bone", "cool", "copper", "flag", "gray", "hot", "hsv",
    "inferno", "jet", "magma", "pink", "plasma", "prism", "spring", 
    "summer", "viridis", "winter",
    
    # GIST colormaps
    "gist_earth", "gist_heat", "gist_ncar", "gist_rainbow", "gist_stern",
    
    # ColorBrewer diverging
    "rBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", 
    "Spectral",
    
    # ColorBrewer sequential
    "Blues", "BuGn", "BuPu", "GnBu", "Greens", "Greys", "Oranges", "OrRd",
    "PuBu", "PuBuGn", "PuRd", "Purples", "RdPu", "Reds", "YlGn", "YlGnBu",
    "YlOrBr", "YlOrRd",
    
    # Miscellaneous
    "afmhot", "brg", "bwr", "coolwarm", "CMRmap", "cubehelix", "gnuplot",
    "gnuplot2", "ocean", "rainbow", "seismic", "nipy_spectral", "terrain"
)

end # module
