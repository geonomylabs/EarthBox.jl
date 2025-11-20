module TestTools

function run_test(test_name::String, TESTS::Dict{String, Function})
    if haskey(TESTS, test_name)
        TESTS[test_name]()
    else
        println("Test '$test_name' not found. Available tests:")
        for name in keys(TESTS)
            println("  - $name")
        end
    end
end


function run_tests(test_names::Vector{String}, TESTS::Dict{String, Function})
    for name in test_names
        run_test(name, TESTS)
    end
end


function run_all_tests(TESTS::Dict{String, Function})
    for test in values(TESTS)
        test()
        println()
    end
end


end