module CaseBuilder

import EarthBox.ParameterRegistry: get_eb_parameters
import ..CaseTypes: CaseCollectionType, CaseParameter

const _FixedScalar = Union{Float64, Integer}
const _FIXED_TRIPLE = Tuple{String, _FixedScalar, String}

function _scalar_for_case_parameter(v::_FixedScalar)
    v isa Bool && return v
    v isa Integer && return Int64(v)
    return v::Float64
end

""" 
    define_case_group!(
        case_inputs::CaseCollectionType;
        case_id_ini::Int,
        parameter_name::String,
        values::Vector{Float64},
        units::String,
        fixed::Union{Nothing, AbstractVector{<:Tuple{String, Union{Float64, Integer}, String}}}=nothing,
        fixed_parameter_name::Union{String, Nothing}=nothing,
        fixed_value::Union{Float64, Integer, Nothing}=nothing,
        fixed_units::Union{String, Nothing}=nothing,
        fixed_parameter_name2::Union{String, Nothing}=nothing,
        fixed_value2::Union{Float64, Integer, Nothing}=nothing,
        fixed_units2::Union{String, Nothing}=nothing,
        fixed_parameter_name3::Union{String, Nothing}=nothing,
        fixed_value3::Union{Float64, Integer, Nothing}=nothing,
        fixed_units3::Union{String, Nothing}=nothing,
        fixed_parameter_name4::Union{String, Nothing}=nothing,
        fixed_value4::Union{Float64, Integer, Nothing}=nothing,
        fixed_units4::Union{String, Nothing}=nothing,
        fixed_parameter_name5::Union{String, Nothing}=nothing,
        fixed_value5::Union{Float64, Integer, Nothing}=nothing,
        fixed_units5::Union{String, Nothing}=nothing,
        fixed_parameter_name6::Union{String, Nothing}=nothing,
        fixed_value6::Union{Float64, Integer, Nothing}=nothing,
        fixed_units6::Union{String, Nothing}=nothing,
        fixed_parameter_name7::Union{String, Nothing}=nothing,
        fixed_value7::Union{Float64, Integer, Nothing}=nothing,
        fixed_units7::Union{String, Nothing}=nothing,
        fixed_parameter_name8::Union{String, Nothing}=nothing,
        fixed_value8::Union{Float64, Integer, Nothing}=nothing,
        fixed_units8::Union{String, Nothing}=nothing,
    )::Int

Define case inputs for a given target parameter name and a list of values. 
Optional fixed keys (up to eight) can be set to fixed values for all cases.

Use either the numbered keyword triples (`fixed_parameter_name`, `fixed_value`,
`fixed_units`, … through `fixed_parameter_name8` / `fixed_value8` / `fixed_units8`)
or the single keyword `fixed` with a vector of `(name, value, units)` tuples, not
both.

# Arguments
- `case_inputs::CaseCollectionType`
    - Dictionary containing the case inputs for each case.
- `case_id_ini::Int`
    - Initial case ID to start building cases from.
- `parameter_name::String`
    - Parameter name in the case inputs dictionary to modify for each case
       using values from the `values` vector.
- `values::Vector{Float64}`
    - List of values to assign to the parameter name for each case.
- `units::String`
    - Units to assign to the parameter name for each case.
- `fixed::Union{Nothing, AbstractVector{<:Tuple{String, Union{Float64, Integer}, String}}}`
    - Optional list of at most eight `(parameter_name, value, units)` triples. When
      not `nothing`, none of the numbered `fixed_parameter_name` / `fixed_value` /
      `fixed_units` keywords may be set. Each `value` may be a `Float64` or any
      `Integer` (for example `Int` or `Int64`).
- `fixed_parameter_name`, `fixed_parameter_name2`, …, `fixed_parameter_name8`
    (`Union{String, Nothing}`)
    - Parameter name keys to set to fixed values for all cases (legacy API).
- `fixed_value`, `fixed_value2`, …, `fixed_value8` (`Union{Float64, Integer, Nothing}`)
    - Fixed values for the corresponding fixed parameter names.
- `fixed_units`, `fixed_units2`, …, `fixed_units8` (`Union{String, Nothing}`)
    - Units for the corresponding fixed parameter names.

# Returns
- `case_id::Int`: The case ID of the last case built

"""
function define_case_group!(
    case_inputs::CaseCollectionType;
    case_id_ini::Int,
    parameter_name::String,
    values::Vector{Float64},
    units::String,
    fixed::Union{Nothing, AbstractVector{<:_FIXED_TRIPLE}}=nothing,
    fixed_parameter_name::Union{String, Nothing}=nothing,
    fixed_value::Union{_FixedScalar, Nothing}=nothing,
    fixed_units::Union{String, Nothing}=nothing,
    fixed_parameter_name2::Union{String, Nothing}=nothing,
    fixed_value2::Union{_FixedScalar, Nothing}=nothing,
    fixed_units2::Union{String, Nothing}=nothing,
    fixed_parameter_name3::Union{String, Nothing}=nothing,
    fixed_value3::Union{_FixedScalar, Nothing}=nothing,
    fixed_units3::Union{String, Nothing}=nothing,
    fixed_parameter_name4::Union{String, Nothing}=nothing,
    fixed_value4::Union{_FixedScalar, Nothing}=nothing,
    fixed_units4::Union{String, Nothing}=nothing,
    fixed_parameter_name5::Union{String, Nothing}=nothing,
    fixed_value5::Union{_FixedScalar, Nothing}=nothing,
    fixed_units5::Union{String, Nothing}=nothing,
    fixed_parameter_name6::Union{String, Nothing}=nothing,
    fixed_value6::Union{_FixedScalar, Nothing}=nothing,
    fixed_units6::Union{String, Nothing}=nothing,
    fixed_parameter_name7::Union{String, Nothing}=nothing,
    fixed_value7::Union{_FixedScalar, Nothing}=nothing,
    fixed_units7::Union{String, Nothing}=nothing,
    fixed_parameter_name8::Union{String, Nothing}=nothing,
    fixed_value8::Union{_FixedScalar, Nothing}=nothing,
    fixed_units8::Union{String, Nothing}=nothing,
)::Int
    keys = get_eb_parameters()
    valid_keys = [getfield(keys, f).name for f in fieldnames(typeof(get_eb_parameters()))]

    if !(parameter_name in valid_keys)
        throw(ArgumentError(
            "Invalid target key: $parameter_name. Expected one of $valid_keys"
        ))
    end

    legacy_triples = (
        (fixed_parameter_name, fixed_value, fixed_units),
        (fixed_parameter_name2, fixed_value2, fixed_units2),
        (fixed_parameter_name3, fixed_value3, fixed_units3),
        (fixed_parameter_name4, fixed_value4, fixed_units4),
        (fixed_parameter_name5, fixed_value5, fixed_units5),
        (fixed_parameter_name6, fixed_value6, fixed_units6),
        (fixed_parameter_name7, fixed_value7, fixed_units7),
        (fixed_parameter_name8, fixed_value8, fixed_units8),
    )
    legacy_any = any(t -> t[1] !== nothing, legacy_triples)
    fixed_kw_used = fixed !== nothing

    if fixed_kw_used && legacy_any
        throw(ArgumentError(
            "Pass either `fixed` or the numbered fixed_parameter_name / fixed_value / " *
            "fixed_units keywords, not both."
        ))
    end

    if fixed_kw_used
        fixed_triples = collect(something(fixed, _FIXED_TRIPLE[]))
    else
        fixed_triples = _FIXED_TRIPLE[]
        for (fpn, fv, fu) in legacy_triples
            if fpn !== nothing
                if fv === nothing
                    throw(ArgumentError(
                        "Fixed value cannot be nothing if fixed key is set ($fpn)."
                    ))
                end
                if fu === nothing
                    throw(ArgumentError(
                        "Fixed units cannot be nothing if fixed key is set ($fpn)."
                    ))
                end
                push!(fixed_triples, (fpn, fv, fu))
            end
        end
    end

    if length(fixed_triples) > 8
        throw(ArgumentError("At most eight fixed parameters are allowed; got $(length(fixed_triples))."))
    end

    seen_fixed = Set{String}()
    for (fpn, fv, fu) in fixed_triples
        if !(fpn in valid_keys)
            throw(ArgumentError(
                "Invalid fixed key: $fpn. Expected one of $valid_keys"
            ))
        end
        if fv === nothing
            throw(ArgumentError("Fixed value cannot be nothing if fixed key is set ($fpn)."))
        end
        if fu === nothing
            throw(ArgumentError("Fixed units cannot be nothing if fixed key is set ($fpn)."))
        end
        if fpn == parameter_name
            throw(ArgumentError(
                "Fixed parameter $fpn cannot match the swept parameter_name."
            ))
        end
        if fpn in seen_fixed
            throw(ArgumentError("Duplicate fixed parameter name: $fpn"))
        end
        push!(seen_fixed, fpn)
    end

    if isempty(values)
        throw(ArgumentError("Values list cannot be empty."))
    end

    icount = 0
    for value in values
        case_id = case_id_ini + icount
        case = "case$case_id"
        case_inputs[case][parameter_name] = CaseParameter(value, units)
        for (fpn, fv, fu) in fixed_triples
            case_inputs[case][fpn] = CaseParameter(_scalar_for_case_parameter(fv), fu)
        end
        icount += 1
    end

    return case_id_ini + icount - 1
end

end # module
