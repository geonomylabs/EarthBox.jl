module GetProperty

using EarthBox.ModelDataContainer: ModelData

function get_radiogenic_heat_production(model::ModelData, material_id::Int16)::Float64
    return model.materials.arrays.mat_hr.array[material_id, 1]
end

function get_heat_capacity_reference(model::ModelData, material_id::Int16)::Float64
    return model.materials.arrays.mat_cp.array[material_id, 1]
end

function get_thermal_conductivity_reference(model::ModelData, material_id::Int)::Float64
    return model.materials.arrays.mat_kt.array[material_id, 1]
end

function get_thermal_conductivity_term_a(model::ModelData, material_id::Int16)::Float64
    return model.materials.arrays.mat_kt.array[material_id, 2]
end

function get_standard_density(model::ModelData, material_id::Int16)::Float64
    return model.materials.arrays.mat_rho.array[material_id,1]
end

function get_thermal_expansivity(model::ModelData, material_id::Int16)::Float64
    return model.materials.arrays.mat_rho.array[material_id, 2]
end

function get_compressibility(model::ModelData, material_id::Int16)::Float64
    return model.materials.arrays.mat_rho.array[material_id, 3]
end

end # module 