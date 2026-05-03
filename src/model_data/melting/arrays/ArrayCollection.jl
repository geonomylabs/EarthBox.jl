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
    Arrays()::Arrays

Create a new Arrays collection with default extraction arrays and empty
scratch buffers. Buffers are sized lazily on first use; see [`Buffers`](@ref).

"""
mutable struct Arrays <: AbstractArrayCollection
    extraction::Extraction
    buffers::Buffers
end

function Arrays()::Arrays
    return Arrays(Extraction(), Buffers())
end

end # module