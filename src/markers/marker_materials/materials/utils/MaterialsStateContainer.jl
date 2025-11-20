module MaterialsStateContainer

struct MaterialsState
    state::String
end

MaterialsState() = MaterialsState("Unloaded")

function update(materials_state::MaterialsState, new_state::String)
    if new_state in get_options()
        return MaterialsState(new_state)
    end
    return materials_state
end

function isloaded(materials_state::MaterialsState)::Bool
    check = false
    if materials_state.state == "Loaded"
        check = true
    end
    return check
end

function get_options()::Vector{String}
    return ["Unloaded", "Loaded"]
end

function check_state(materials_state::MaterialsState)
    println("materials state: $(materials_state.state)")
    if materials_state.state == "Unloaded"
        error(
            """
            Materials have not been loaded and, therefore, marker plots 
            cannot be created. Please define material_library_file_path 
            and material_model_file_path or materials_input_dict during 
            ModelPlots2D initialization. See the method initialization() 
            in ModelPlots2D for details.
            """
            )
    end
end

end # module 