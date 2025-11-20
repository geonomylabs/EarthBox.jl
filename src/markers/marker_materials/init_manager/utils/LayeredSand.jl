module LayeredSand

"""
    define_layered_structure_in_sand(
        matid::Int,
        matid_layered_sand_b::Int,
        y_marker::Float64,
        ysize::Float64,
        nsand_layers::Int
    )::Int

Define layered structure in sand.

# Arguments
- `matid`: Material ID of marker
- `matid_layered_sand_b`: Material ID of layered sand B
- `y_marker`: Y-coordinate of marker
- `ysize`: Y-size of model
- `nsand_layers`: Number of sand layers

# Returns
- Material ID of marker
"""
function define_layered_structure_in_sand(
    matid::Int16,
    matid_layered_sand_b::Int16,
    y_marker::Float64,
    ysize::Float64,
    nsand_layers::Int
)::Int16
    m2_check = ceil(y_marker/(ysize/nsand_layers) - 0.5) + 1
    m2_check = m2_check - ceil(m2_check/2 - 0.5)*2
    if m2_check == 0
        matid = matid_layered_sand_b
    end
    return matid
end

end # module LayeredSand 