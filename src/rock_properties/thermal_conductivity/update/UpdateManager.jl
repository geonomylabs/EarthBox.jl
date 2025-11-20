module UpdateManager

include("../types/Liao14.jl")
include("../types/Uniform.jl")
include("../types/TemperatureDependent.jl")
include("../types/Naliboff17.jl")
include("../types/Adiabatic.jl")
include("../types/SekiguchiWaples.jl")

import EarthBox.ModelDataContainer: ModelData
import EarthBox.PrintFuncs: @timeit_memit
import EarthBox: OptionTools
import ..Options: option_names
import ..Options: option_ids
import .Liao14
import .Uniform
import .TemperatureDependent
import .Naliboff17
import .Adiabatic
import .SekiguchiWaples

function update_conductivity!(
    model::ModelData, 
    inside_flags::Vector{Int8}, 
    ::Val{option_names.Liao14}
)::Nothing
    @timeit_memit "Finished updating marker thermal conductivity using Liao14 model" begin
        Liao14.update_conductivity!(model, inside_flags)
    end
end

function update_conductivity!(
    model::ModelData, 
    inside_flags::Vector{Int8}, 
    ::Val{option_names.Uniform}
)::Nothing
    @timeit_memit "Finished updating marker thermal conductivity using Uniform model" begin
        Uniform.update_conductivity!(model, inside_flags)
    end
end

function update_conductivity!(
    model::ModelData, 
    inside_flags::Vector{Int8}, 
    ::Val{option_names.TemperatureDependent}
)::Nothing
    @timeit_memit "Finished updating marker thermal conductivity using TemperatureDependent model" begin
        TemperatureDependent.update_conductivity!(model, inside_flags)
    end
end

function update_conductivity!(
    model::ModelData, 
    inside_flags::Vector{Int8}, 
    ::Val{option_names.Naliboff17}
)::Nothing
    @timeit_memit "Finished updating marker thermal conductivity using Naliboff17 model" begin
        Naliboff17.update_conductivity!(model, inside_flags)
    end
end

function update_conductivity!(
    model::ModelData, 
    inside_flags::Vector{Int8}, 
    ::Val{option_names.Adiabatic}
)::Nothing
    @timeit_memit "Finished updating marker thermal conductivity using Adiabatic model" begin
        Adiabatic.update_conductivity!(model, inside_flags)
    end
end

function update_conductivity!(
    model::ModelData, 
    inside_flags::Vector{Int8}, 
    ::Val{option_names.SekiguchiWaples}
)::Nothing
    @timeit_memit "Finished updating marker thermal conductivity using SekiguchiWaples model" begin
        SekiguchiWaples.update_conductivity!(model, inside_flags)
    end
end

end # module
