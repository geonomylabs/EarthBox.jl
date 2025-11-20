module MeltCheck

""" Check if melt is available for extraction.

# Returns
- check::Bool
    - True if melt is available for extraction.
"""
function melt_available(
    meltfrac::Float64,
    extractable_meltfrac::Float64
)::Bool
    check = false
    if meltfrac > 0 && extractable_meltfrac > 0
        check = true
    end
    return check
end

end # module 