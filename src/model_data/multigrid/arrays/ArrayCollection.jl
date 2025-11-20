module ArrayCollection

include("groups/__init__.jl")

import EarthBox.EarthBoxDtypes: AbstractArrayCollection
import .SmoothingGroup: Smoothing
import .MeanResidualsGroup: MeanResiduals

mutable struct Arrays <: AbstractArrayCollection
    smoothing::Smoothing
    mean_residuals::MeanResiduals
end

function Arrays(
    nvcycles::Int64, 
    smoothing_iterations_on_finest_level::Int64,
    max_levels::Int64 
)::Arrays
    return Arrays(
        Smoothing(smoothing_iterations_on_finest_level, max_levels),
        MeanResiduals(nvcycles)
    )
end

end # module 