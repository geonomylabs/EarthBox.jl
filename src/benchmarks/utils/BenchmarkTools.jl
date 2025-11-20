module BenchmarkTools

import ..Profiles
import ..BenchmarksStruct: Benchmarks

""" Extract values from numerical models along x-direction profiles.
"""
function get_xprofile_from_numerical_model(
    bench::Benchmarks,
    file_base_name::String
)::Tuple{Vector{Float64}, Vector{Float64}, Float64}
    get_xprofile_func_dict = Dict(
        "vel_cmyr" => Profiles.get_xprofile_vy,
        "log_eta_Pas" => Profiles.get_xprofile_scalar,
        "TempC" => Profiles.get_xprofile_scalar
    )

    x_coors, numerical_values, tmyr = get_xprofile_func_dict[file_base_name](
        bench.main_paths["post_proc_input_path"],
        file_base_name,
        bench.itime_step,
        bench.test_info_dict[bench.main_paths["model_name"]]["y_index"]
    )

    return x_coors, numerical_values, tmyr
end

end # module 