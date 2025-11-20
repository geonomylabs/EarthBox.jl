module ParameterCollection

include("groups/__init__.jl")

import EarthBox.EarthBoxDtypes: AbstractParameterCollection
import .VcycleGroup: Vcycle
import .RelaxationGroup: Relaxation
import .RenormalizeViscosityGroup: RenormalizeViscosity

mutable struct Parameters <: AbstractParameterCollection
    vcycle::Vcycle
    relaxation::Relaxation
    renormalize_viscosity::RenormalizeViscosity
end

function Parameters(nvcycles::Int64, levelnum::Int64, max_levels::Int64)::Parameters
    return Parameters(
        Vcycle(nvcycles, levelnum, max_levels),
        Relaxation(),
        RenormalizeViscosity(),
    )
end

end # module 