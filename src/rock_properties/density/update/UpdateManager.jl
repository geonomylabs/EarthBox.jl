module UpdateManager

include("../types/DensityBlankenbach89.jl")
include("../types/DensityLiao14.jl")
include("../types/DensityNaliboff17.jl")

import EarthBox.ModelDataContainer: ModelData
import EarthBox.PrintFuncs: @timeit_memit
import EarthBox: OptionTools
import ..Options: option_names
import ..Options: option_ids
import .DensityBlankenbach89
import .DensityLiao14
import .DensityNaliboff17

function update_density!(
    model::ModelData,
    inside_flags::Vector{Int8},
    ::Val{option_names.Blankenbach89}
)::Nothing
    @timeit_memit "Finished updating marker density using Blankenbach89 model" begin
        DensityBlankenbach89.update_density!(model, inside_flags)
    end
end

function update_density!(
    model::ModelData,
    inside_flags::Vector{Int8},
    ::Val{option_names.Liao14}
)::Nothing
    @timeit_memit "Finished updating marker density using Liao14 model" begin
        DensityLiao14.update_density!(model, inside_flags)
    end
end

function update_density!(
    model::ModelData,
    inside_flags::Vector{Int8},
    ::Val{option_names.Naliboff17}
)::Nothing
    @timeit_memit "Finished updating marker density using Naliboff17 model" begin
        DensityNaliboff17.update_density!(model, inside_flags)
    end
end

end # module
