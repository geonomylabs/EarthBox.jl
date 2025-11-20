module UpdateManager

include("../types/ConstantRhoCp.jl")
include("../types/TemperatureDependentWaples.jl")

import EarthBox.ModelDataContainer: ModelData
import EarthBox.PrintFuncs: @timeit_memit
import ..Options: option_names
import .ConstantRhoCp
import .TemperatureDependentWaples

function update_rhocp!(
    model::ModelData,
    inside_flags::Vector{Int8}, 
    ::Val{option_names.Constant}
)::Nothing
    @timeit_memit "Finished updating marker rhocp using ConstantRhoCp model" begin
        ConstantRhoCp.update_rhocp!(model, inside_flags)
    end
    return nothing
end

function update_rhocp!(
    model::ModelData, 
    inside_flags::Vector{Int8}, 
    ::Val{option_names.TemperatureDependentWaples}
)::Nothing
    @timeit_memit "Finished updating marker rhocp using TemperatureDependentWaples model" begin
        TemperatureDependentWaples.update_rhocp!(model, inside_flags)
    end
    return nothing
end

end # module UpdateManager
