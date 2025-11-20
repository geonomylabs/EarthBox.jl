module LayerMaterialIDsContainer

mutable struct LayerMaterialIDs
    # Maximum number of material ID's
    nid_max::Int
    # Material ID's of sticky air
    matids_sticky::Vector{Int16}
    # Material ID's of sticky water
    matids_sticky_water::Vector{Int16}
    # Material ID's of rock
    matids_rock::Vector{Int16}
    # Material ID's of crust
    matids_crust::Vector{Int16}
    # Material ID's of lithosphere
    matids_lithosphere::Vector{Int16}
    # Material ID's of asthenosphere
    matids_asthenosphere::Vector{Int16}
    # Material ID's of sediments
    matids_sediments::Vector{Int16}
end

# Constructor with default values
function LayerMaterialIDs()
    # Maximum number of material ID's
    nid_max = 20
    # Initialize arrays with -1
    matids_sticky = fill(Int16(-1), nid_max)
    matids_sticky_water = fill(Int16(-1), nid_max)
    matids_rock = fill(Int16(-1), nid_max)
    matids_crust = fill(Int16(-1), nid_max)
    matids_lithosphere = fill(Int16(-1), nid_max)
    matids_asthenosphere = fill(Int16(-1), nid_max)
    matids_sediments = fill(Int16(-1), nid_max)
    
    return LayerMaterialIDs(nid_max, matids_sticky, matids_sticky_water,
        matids_rock, matids_crust, matids_lithosphere, matids_asthenosphere,
        matids_sediments)
end

function set_material_ids!(
    data::LayerMaterialIDs,
    matids_sticky::Vector{Int16},
    matids_sticky_water::Vector{Int16},
    matids_rock::Vector{Int16},
    matids_crust::Vector{Int16},
    matids_lithosphere::Vector{Int16},
    matids_asthenosphere::Vector{Int16},
    matids_sediments::Union{Vector{Int16}, Nothing}=nothing
)::Nothing
    for (i, matid) in enumerate(matids_sticky)
        if i <= data.nid_max
            data.matids_sticky[i] = matid
        end
    end
    
    for (i, matid) in enumerate(matids_sticky_water)
        if i <= data.nid_max
            data.matids_sticky_water[i] = matid
        end
    end
    
    for (i, matid) in enumerate(matids_rock)
        if i <= data.nid_max
            data.matids_rock[i] = matid
        end
    end
    
    for (i, matid) in enumerate(matids_crust)
        if i <= data.nid_max
            data.matids_crust[i] = matid
        end
    end
    
    for (i, matid) in enumerate(matids_lithosphere)
        if i <= data.nid_max
            data.matids_lithosphere[i] = matid
        end
    end
    
    for (i, matid) in enumerate(matids_asthenosphere)
        if i <= data.nid_max
            data.matids_asthenosphere[i] = matid
        end
    end
    
    if isnothing(matids_sediments)
        data.matids_sediments = fill(-1, data.nid_max)
    else
        for (i, matid) in enumerate(matids_sediments)
            if i <= data.nid_max
                data.matids_sediments[i] = matid
            end
        end
    end
    return nothing
end

end # module 