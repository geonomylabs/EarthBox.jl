module TestResults

using Printf: @sprintf

""" Check test results by comparing numerical and analytical solutions.
"""
function get_test_results_numerical_vs_analytical(
    numerical_solution::Vector{Float64},
    analytical_solution::Vector{Float64},
    relative_error_limit::Float64
)::Tuple{Vector{Union{String,Float64}},String}
    difference = abs.(numerical_solution .- analytical_solution)
    relative_error = calc_relative_error_1d(
        analytical_solution,
        numerical_solution
    )
    max_relative_error = maximum(relative_error)

    print_benchmark_error(
        difference, max_relative_error, relative_error_limit
    )

    if max_relative_error <= relative_error_limit
        result = [
            "Success",
            max_relative_error * 100.0,
            relative_error_limit * 100.0
        ]
        result_msg = string(
            "Test Successful: max relative error is not greater than limit: ",
            "$(@sprintf("%.3f", max_relative_error*100.0))% <= ",
            "$(@sprintf("%.3f", relative_error_limit*100.0))%"
        )
    else
        result = [
            "Failure",
            max_relative_error * 100.0,
            relative_error_limit * 100.0
        ]
        result_msg = string(
            "!! Test Failed !! : max relative error is greater than limit: ",
            "$(@sprintf("%.3f", max_relative_error*100.0))% > ",
            "$(@sprintf("%.3f", relative_error_limit*100.0))%"
        )
    end
    return result, result_msg
end

""" Calculate the relative error between 1D vectors.
"""
function calc_relative_error_1d(
    analytical_vector::Vector{Float64},
    numerical_vector::Vector{Float64}
)::Vector{Float64}
    nvalues = length(analytical_vector)
    relative_error = zeros(Float64, nvalues)
    for i in 1:nvalues
        value_analytical = analytical_vector[i]
        value_numerical = numerical_vector[i]
        difference = abs(value_numerical - value_analytical)
        if abs(value_analytical) > 0
            relative_error[i] = difference / abs(value_analytical)
        end
    end
    return relative_error
end

""" Print information on benchmark error.
"""
function print_benchmark_error(
    difference::Vector{Float64},
    max_relative_error::Float64,
    relative_error_limit::Float64
)::Nothing
    println("")
    println("Max difference (num-ana): $(@sprintf("%.5f", maximum(difference)))")
    println("Max relative error: $(@sprintf("%.3f", max_relative_error*100.0)) %")
    println("Max relative error limit: $(@sprintf("%.3f", relative_error_limit*100.0)) %")
    return nothing
end

""" Print test result.
"""
function print_test_result(result::Any)::Nothing
    println("")
    println(result)
    return nothing
end

end # module 