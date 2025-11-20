module Collections

mutable struct TransportCollections
    use_collections::Bool
    topo_collection::Dict{Int, Vector{Float64}}
    water_depth_collection::Dict{Int, Vector{Float64}}
    divides_collection::Dict{Int, Vector{Float64}}
    basement_collection::Dict{Int, Vector{Float64}}
end

function TransportCollections(use_collections::Bool=true)
    return TransportCollections(
        use_collections,
        Dict{Int, Vector{Float64}}(),
        Dict{Int, Vector{Float64}}(),
        Dict{Int, Vector{Float64}}(),
        Dict{Int, Vector{Float64}}()
    )
end

function update_collections!(
    collections::TransportCollections, 
    istep::Int, 
    topo_gridy::Vector{Float64}, 
    drainage_divides_x::Vector{Float64},
    water_depth_x::Vector{Float64}
)
    if collections.use_collections
        collections.topo_collection[istep] = copy(topo_gridy)
        collections.divides_collection[istep] = copy(drainage_divides_x)
        collections.water_depth_collection[istep] = copy(water_depth_x)
        collections.basement_collection[istep] = copy(topo_gridy)
    end
end

end