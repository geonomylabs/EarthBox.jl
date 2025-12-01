module BenchmarksManager

include("utils/BenchmarksStruct.jl")
include("utils/Reader.jl")
include("utils/Profiles.jl")
include("utils/BenchmarkTools.jl")
include("utils/BenchmarkFuncs.jl")
include("utils/TestResults.jl")
include("options/Options.jl")
include("post_processing/PostProcessingManager.jl")

using YAML
import Dates
import DataStructures: OrderedDict
import EarthBox
import EarthBox: PrintFuncs
import EarthBox.EarthBoxPaths: EarthBoxPathsState
import EarthBox.MaterialLibraryCollection: MaterialLibrary
import EarthBox.OptionTools: get_id
import EarthBox: ModelManager
import EarthBox.PrintFuncs: PRINT_SETTINGS
import .Options: OptionState
import .Options: get_options
import .Options: option_ids
import .Options: option_names
import .BenchmarksStruct: Benchmarks

export get_options, option_names, run_benchmarks, run_benchmark
export couette_flow_viscous_heating

BM_OPTIONS = get_options()

function get_benchmark_descriptions_string()::String
    return """
    # Benchmark Descriptions:
    $(join(["`:$(name)` \n - $(BM_OPTIONS[option_ids[name]].description)" for name in option_names],"\n\n"))
    """
end

function get_mumps_args_string()::String
    return """
    # Additional Keyword Arguments:
    - `use_mumps::Bool = false`:
        - Set to true to use the MUMPS solver. Default is false.
    - `nprocs::Int = 1`:
        - Number of processors to use if using the MUMPS solver. Default is 1.
    """
end

function get_keyword_args_string()::String
    return """
    # Optional Keyword Arguments:
    - `run_model::Bool`:
        - Set to true to run the model. Default is true.
    - `run_post_processing::Bool`:
        - Set to true to run post processing. Default is true.
    - `base_path::String`:
        - Base path for the output directory. EarthBox creates a time-stamped subdirectory in this path. 
           Default is the present working directory (pwd()).
    - `old_date_stamp::String`:
        - Old date stamp for the output directory. If equal to nothing, a new date stamp will be generated. 
           If a old time stamp is provided, the output directory will use the old time stamp.
           This is useful if you want to re-run post processing on a previous run. The format of the time stamp is
           YYYY-MM-DD_HH-MM-SS. For example: "2024-04-20_12-53-59". Default is nothing.
    - `make_backup::Bool`:
        - Set to true to make backup model files with each output step. This is useful for testing 
           model restart functionality. Default is false.
    - `restart_from_backup::Bool`:
        - Set to true to restart from a backup file. Default is false.
    
    # Returns
    - Nothing
    """
end

"""
    run_benchmark(benchmark_name::Symbol; kwargs...)::Nothing

Run the benchmark specified by the benchmark_name.

$(get_keyword_args_string())

$(get_mumps_args_string())

# Example:

To run the `:couette_flow_viscous_heating` benchmark with output being sent to the present working 
directory:

```julia
using EarthBox
BenchmarksManager.run_benchmark(:couette_flow_viscous_heating);
```

$(get_benchmark_descriptions_string())
"""
function run_benchmark(benchmark_name::Symbol; kwargs...)::Tuple{String, Float64, Float64}
    valid_options = [v for v in option_names]
    if !(benchmark_name in valid_options)
        error("Invalid benchmark name: $benchmark_name. Valid benchmark names are: $(valid_options)")
    end
    test_dict = OrderedDict(benchmark_name => true)
    bench = run_benchmarks(test_dict; mumps_solver_dict=get_mumps_solver_dict(kwargs...), kwargs...)

    results_dict = bench.results_dict
    test_results_dict = results_dict[String(benchmark_name)]["summary_for_each_time_step"]
    last_key = collect(keys(test_results_dict))[end]
    result_dict = test_results_dict[last_key]

    status_str = result_dict["result"] # "Success" or "Failure"
    max_relative_error_percentage = result_dict["max relative error percentage"]
    relative_error_limit_percentage = result_dict["relative error limit %"]

    return status_str, max_relative_error_percentage, relative_error_limit_percentage
end

"""
    run_benchmarks(
        test_dict::OrderedDict{Symbol, Bool};
        mumps_solver_dict::Union{Dict{Symbol, Vector{Any}}, Nothing} = nothing,
        kwargs...
    )::Benchmarks

Run the benchmarks specified in the test_dict.

# Arguments
- `test_dict::OrderedDict{Symbol, Bool}`:
    - Dictionary of benchmarks to run. The key is the benchmark name and the 
       value is a boolean indicating if the benchmark should be run.
- `mumps_solver_dict::Union{Dict{Symbol, Vector{Any}}, Nothing} = nothing`:
    - Dictionary of benchmarks to run with MUMPS solver. The key is the benchmark 
       name and the value is a vector of two elements. The first element is a boolean 
       indicating if the benchmark should be run with MUMPS solver and the second 
       element is the number of processors to use.

$(get_keyword_args_string())

# Returns
- `bench::Benchmarks`:
    - Benchmarks object containing the results of the benchmarks.

"""
function run_benchmarks(
    test_dict::OrderedDict{Symbol, Bool};
    mumps_solver_dict::Union{Dict{Symbol, Vector{Any}}, Nothing} = nothing,
    kwargs...
)::Benchmarks
    PRINT_SETTINGS.print_performance = false
    PRINT_SETTINGS.print_info = false

    run_model = get(kwargs, :run_model, true) 
    run_post_processing = get(kwargs, :run_post_processing, true)
    base_path = get(kwargs, :base_path, pwd())
    old_date_stamp = get(kwargs, :old_date_stamp, nothing)
    make_backup = get(kwargs, :make_backup, false)
    restart_from_backup = get(kwargs, :restart_from_backup, false)
    output_dir_root = get_output_dir_root(base_path, old_date_stamp=old_date_stamp)
    bench = Benchmarks(
        run_model           = run_model,
        run_post_processing = run_post_processing,
        output_dir_root     = output_dir_root
    )
    run_tests(
        bench, test_dict, mumps_solver_dict;
        make_backup         = make_backup, 
        restart_from_backup = restart_from_backup,
    )
    return bench
end

function run_tests(
    bench::Benchmarks,
    test_dict::OrderedDict{Symbol, Bool},
    mumps_solver_dict::Union{Dict{Symbol, Vector{Any}}, Nothing} = nothing;
    make_backup::Bool = false,
    restart_from_backup::Bool = false
)::Nothing
    if isnothing(test_dict)
        error("test_dict is not defined.")
    end
    for (test_name, isactive) in test_dict
        if !isa(test_name, Symbol)
            error("test_dict key $test_name is not a Symbol. Correct inputs")
        end
        if isactive
            println("\n=============================================================")
            println("Working on benchmark test: $test_name")
            println("=============================================================\n")
            option = get_option_state(String(test_name))
            PrintFuncs.print_option_info(String(test_name), option.description)
            bench.main_paths["model_name"] = String(test_name)
            if bench.run_model
                if !isnothing(mumps_solver_dict) && haskey(mumps_solver_dict, test_name)
                    use_mumps = Bool(mumps_solver_dict[test_name][1])
                    nprocs = mumps_solver_dict[test_name][2]
                else
                    use_mumps = false
                    nprocs = 1
                end
                run_test_model(
                    bench, String(test_name), use_mumps, nprocs;
                    make_backup = make_backup,
                    restart_from_backup = restart_from_backup
                )
            end
            if bench.run_post_processing
                run_post_test_processing(bench, String(test_name))
            end
        end
    end
    return nothing
end

function run_test_model(
    bench::Benchmarks,
    test_name::String,
    use_mumps::Bool,
    nprocs::Int;
    make_backup::Bool = false,
    restart_from_backup::Bool = false
)::Nothing
    model_dir = joinpath(bench.main_paths["models_path"], test_name)
    ebpaths = EarthBoxPathsState().key_names

    eb = EarthBox.EarthBoxState(
        restart_from_backup = restart_from_backup,
        paths = Dict(
            ebpaths.model_input_file => joinpath(model_dir, "model.yml"),
            ebpaths.materials_input_file => joinpath(model_dir, "materials.yml"),
            ebpaths.materials_library_file => bench.main_paths["materials_library_file"],
            ebpaths.output_dir => joinpath(
                bench.main_paths["output_dir_root"], test_name * "_output"
            )
        ),
        use_mumps = use_mumps,
        nprocs = nprocs,
        use_internal_mumps = true
    )

    ModelManager.initialize_model!(eb.model_manager)
    EarthBox.run_time_steps(eb, make_backup = make_backup)
    return nothing
end

function run_post_test_processing(bench::Benchmarks, test_name::String)
    (
        bench.main_paths["post_proc_input_path"],
        bench.main_paths["post_proc_output_path"]
    ) = get_post_processing_paths(bench, test_name)

    bench.model_dict = get_model_dict(bench.main_paths, test_name)
    bench.materials_dict = get_materials_dict(bench.main_paths, test_name)

    test_results_dict = execute_post_processing(bench)

    bench.results_dict[test_name] = Dict(
        "summary_for_each_time_step" => test_results_dict,
        "benchmark_description" => replace(
            get_option_state(test_name).description,
            "\n" => ""
        )
    )

    YAML.write_file(bench.main_paths["result_file_path"], bench.results_dict)
    test_assertion(test_results_dict)

    print_output_location_message(bench.main_paths)
end

function get_post_processing_paths(
    bench::Benchmarks,
    test_name::String
)::Tuple{String, String}
    post_proc_input_path = joinpath(
        bench.main_paths["output_dir_root"],
        test_name * "_output"
    )
    post_proc_output_path = joinpath(post_proc_input_path, "benchmarks")
    return post_proc_input_path, post_proc_output_path
end

function execute_post_processing(
    bench::Benchmarks
)::Dict{String, Dict{String, Any}}
    test_name = bench.main_paths["model_name"]
    input_path = bench.main_paths["post_proc_input_path"]
    output_path = bench.main_paths["post_proc_output_path"]

    option_active = get_option_state(test_name)
    time_steps = option_active.time_steps
    y_index = option_active.y_index

    test_info_dict = Dict(
        test_name => Dict(
            "time_steps" => time_steps,
            "y_index" => y_index
        )
    )

    bench.test_info_dict = test_info_dict

    println(">> Working on post processing...")
    println("    ++ input_path: ", input_path)
    println("    ++ output_path: ", output_path)

    itime_step_max = maximum(time_steps)
    test_results_dict = Dict{String, Dict{String, Any}}()

    for itime_step in time_steps
        println(">> Working on time step $itime_step")
        bench.itime_step = itime_step
        result, result_msg = PostProcessingManager.post_proc_func(
            bench, Val(Symbol(test_name)))

        # Only record all time steps for non-steady benchmarks
        if test_name in [
            String(option_names.channel_flow_non_steady_temperature),
            String(option_names.solid_body_rotation),
            String(option_names.elastic_slab)
        ]
            TestResults.print_test_result(result_msg)
            test_results_dict["time_step_$itime_step"] = Dict(
                "result" => result[1],
                "max relative error percentage" => round(result[2], digits=4),
                "relative error limit %" => round(result[3], digits=4),
                "message" => result_msg
            )
        else
            if itime_step == itime_step_max
                TestResults.print_test_result(result_msg)
                test_results_dict["time_step_$itime_step"] = Dict(
                    "result" => result[1],
                    "max relative error percentage" => round(result[2], digits=4),
                    "relative error limit %" => round(result[3], digits=4),
                    "message" => result_msg
                )
            end
        end
    end

    return test_results_dict
end

function test_assertion(test_results_dict::Dict{String, Dict{String, Any}})
    for (_, result_dict) in test_results_dict
        result = result_dict["result"]
        msg = result_dict["message"]
        @assert result == "Success" msg
    end
end

function get_model_dict(
    main_paths::Dict{String, String},
    test_name::String
)::Dict{String, Any}
    model_path = joinpath(main_paths["models_path"], test_name, "model.yml")
    return YAML.load_file(model_path)
end

function get_materials_dict(
    main_paths::Dict{String, String},
    test_name::String
)::Dict{Any, Any}
    materials_path = joinpath(
        main_paths["models_path"], test_name, "materials.yml"
    )
    return YAML.load_file(materials_path)
end

function print_output_location_message(main_paths::Dict{String, String})
    post_proc_input_path = main_paths["post_proc_input_path"]
    println("\nLook at the following path for benchmark plots: $post_proc_input_path \n")
end

function get_option_state(test_name::String)::OptionState
    options = get_options()
    option_id = get_id(options, test_name)
    return get_options()[option_id]
end

"""
Get output directory root.

Inputs
------
base_path : String
    Base path for output directory.

old_date_stamp : Union{String, Nothing}
    Date stamp for output directory. If equal to nothing, a new date stamp will be generated.

Returns
-------
output_dir_root : String
    Output directory root.
"""
function get_output_dir_root(
    base_path::String;
    old_date_stamp::Union{String, Nothing}=nothing
)::String
    if isnothing(old_date_stamp)
        output_dir_root = joinpath(
            base_path,
            "earthbox_benchmark_results_" * get_date_string()
        )
    else
        output_dir_root = joinpath(
            base_path,
            "earthbox_benchmark_results_" * old_date_stamp
        )
    end
    return output_dir_root
end

"""
Get date string.
"""
function get_date_string()::String
    return Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
end

struct ValidInputNames
    run_model::Symbol
    run_post_processing::Symbol
    base_path::Symbol
    old_date_stamp::Symbol
    make_backup::Symbol
    restart_from_backup::Symbol
    nprocs::Symbol
    use_mumps::Symbol
end

function validate_input_names(kwargs::Dict{Symbol, Any})::Nothing
    valid_names = fieldnames(ValidInputNames)
    for (key, _) in kwargs
        if !(key in valid_names)
            error("Invalid input parameter name: $key. Valid names are: $valid_names")
        end
    end
end

function get_mumps_solver_dict(;kwargs...)::Union{Dict{Symbol, Vector{Any}}, Nothing}
    use_mumps = get(kwargs, :use_mumps, false)
    nprocs = get(kwargs, :nprocs, 1)
    if use_mumps
        mumps_solver_dict = Dict(
            option_names.flexure_triangular_hole => [true, nprocs]
        )
    else
        mumps_solver_dict = nothing
    end
    return mumps_solver_dict
end

end # module