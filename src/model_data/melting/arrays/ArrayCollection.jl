module ArrayCollection

include("groups/__init__.jl")

import EarthBox.EarthBoxDtypes: AbstractArrayCollection
import .ExtractionGroup: Extraction

"""
    Arrays <: AbstractArrayCollection

Data structure containing array groups for melt extraction tracking.

# Fields
- `extraction::`[`Extraction`](@ref): Arrays for drainage basin tracking

# Constructor
    Arrays()::Arrays

Create a new Arrays collection with default extraction arrays.

"""
mutable struct Arrays <: AbstractArrayCollection
    extraction::Extraction
end

function Arrays()::Arrays
    return Arrays(Extraction())
end

end # module 