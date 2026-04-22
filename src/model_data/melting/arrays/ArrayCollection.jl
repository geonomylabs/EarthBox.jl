module ArrayCollection

include("groups/__init__.jl")

import EarthBox.EarthBoxDtypes: AbstractArrayCollection
import .ExtractionGroup: Extraction
import .BuffersGroup: Buffers

"""
    Arrays <: AbstractArrayCollection

Data structure containing array groups for melt extraction tracking.

# Fields
- `extraction::`[`Extraction`](@ref): Arrays for drainage basin tracking
- `buffers::`[`Buffers`](@ref): Pre-allocated marker-sized scratch buffers

# Constructor
    Arrays(marknum::Int)::Arrays

Create a new Arrays collection with default extraction arrays and
marker-sized scratch buffers.

"""
mutable struct Arrays <: AbstractArrayCollection
    extraction::Extraction
    buffers::Buffers
end

function Arrays(marknum::Int)::Arrays
    return Arrays(Extraction(), Buffers(marknum))
end

end # module