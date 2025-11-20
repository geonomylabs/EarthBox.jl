module ViscoelasticExtension

import EarthBox.EarthBoxDtypes
import ...BenchmarksStruct: Benchmarks
import ..ViscoelasticLithosphericDeformation: execute_postprocessing_steps

function compare_numerical_to_analytical(
    bench::Benchmarks
)::Tuple{Vector{Union{String, Float64}}, String}
    result, result_msg = execute_postprocessing_steps(
        bench, "Viscoelastic Extension Benchmark", (-2.0, 10.0))
    return result, result_msg
end

end # module 