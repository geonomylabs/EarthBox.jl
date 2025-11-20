module TupleTools

function merge_named_tuples(tuples...)
    result = tuples[1]
    ntuples = length(tuples)
    for i in 2:ntuples
        result = merge(result, tuples[i])
    end
    return result
end

end # module