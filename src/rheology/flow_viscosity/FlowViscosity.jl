module FlowViscosity

include("options/Options.jl")
include("update/UpdateManager.jl")

import EarthBox.PrintFuncs: print_info
import EarthBox.ModelDataContainer: ModelData
import EarthBox.OptionTools: OptionState
import EarthBox.OptionTools: get_option_symbol_from_id
import .Options: option_names
import .Options: option_ids
import .Options: get_options
import .UpdateManager

""" Update marker flow viscosity.

The flow viscosity array ``marker_eta_flow`` is updated for all markers
in the model using the flow law associated with the markers material id.

The maximum integer ID value of the material flow type array
``mat_flow_type`` is used to determine the flow law function called by
the model. See the get_options() method for a description of each flow
law associated with each integer ID.

Although the maximum flow law ID controls the flow law function,
multiple flow law ID's may be used in the function. For example, if the
maximum ID is 1 or 2 the composite rheology flow law function is called
and the isoviscous flow law ID is also used if applicable to a given
marker.

To Do: Composite viscosity ID's should be consolidated to a single ID.

"""
function update_marker_flow_viscosity!(model::ModelData, inside_flags::Vector{Int8})::Nothing
    option_id_dominate = get_dominate_flow_type_id(model)
    option_name = get_option_symbol_from_id(option_ids, option_id_dominate)
    option_name = maintain_backward_compatibility(option_name)
    UpdateManager.update_marker_flow_viscosity!(model, inside_flags, Val(option_name))
    return nothing
end

function get_dominate_option(
    model::ModelData
)::OptionState
    option_id = get_dominate_flow_type_id(model)
    return get_option(option_id)
end

function get_option(
    option_id::Int
)::OptionState
    return get_options()[option_id]
end

function get_dominate_flow_type_id(model::ModelData)::Int
    mat_flow_type = get_material_flow_type_array(model)
    return maximum(mat_flow_type)
end

function get_material_flow_type_array(model::ModelData)::Vector{Int64}
    return model.materials.arrays.mat_flow_type.array
end

function maintain_backward_compatibility(option_name::Symbol)::Symbol
    if option_name === :Composite2
        return :Composite
    else
        return option_name
    end
end

function print_model_flow_type_information(model::ModelData)::Nothing
    print_material_flow_type_information(model)
    return nothing
end

function print_material_flow_type_information(model::ModelData)::Nothing
    print_info("", level=1)
    print_info("Flow law types used in material model:", level=1)

    mat_flow_type = model.materials.arrays.mat_flow_type.array
    flow_itype_max = maximum(mat_flow_type)
    found_types = Int[]
    
    for flow_itype in mat_flow_type
        if !(flow_itype in found_types) && is_applicable_flow_type(flow_itype)
            print_flow_law_info(flow_itype)
            push!(found_types, flow_itype)
            print_info("", level=1)
        end
    end
    
    return nothing
end

function is_applicable_flow_type(flow_itype::Int)::Bool
    return flow_itype != -9999
end

function print_flow_law_info(flow_itype::Int)::Nothing
    name = String(get_option(flow_itype).option_name)
    description = get_option(flow_itype).description
    print_info("flow_itype = $flow_itype: $name", level=2)
    print_info("-"^79, level=2)
    print_info(description, level=2)
    return nothing
end

end # module FlowViscosity

