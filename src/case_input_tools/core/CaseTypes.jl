module CaseTypes

import DataStructures: OrderedDict

"""
    CaseParameter

Type for storing a model case parameter value and its units.

EarthBox will convert the case parameter value to standard units when the 
parameters are loaded into the model.

"""
struct CaseParameter
    value::Union{Float64, Int64, Bool, String, Symbol}
    units::String
end

"""
    CaseType

Dictionary type for storing model case inputs where the key is the earthbox 
parameter name and the value is a CaseParameter type.
"""
const CaseType = Dict{String, CaseParameter}

"""
    CaseCollectionType

Ordered dictionary type for storing a collection of model cases where the key 
is the case name (`case0`, `case1`, etc) and the value is a CaseType dictionary.
"""
const CaseCollectionType = OrderedDict{String, CaseType}

end # module