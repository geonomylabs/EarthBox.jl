module BenchmarksStruct

import Dates
import EarthBox.EarthBoxDtypes: MaterialsDictType
import EarthBox.EarthBox.EarthBoxDtypes: ModelInputDictType
import EarthBox.MaterialLibraryCollection: MaterialLibrary
import EarthBox.Markers.MarkerMaterials.MaterialsContainer: MaterialReader

mutable struct Benchmarks
    run_model::Bool
    run_post_processing::Bool
    test_dir::String
    project_path::String
    output_dir_root::Union{String, Nothing}
    materials_library_file::String
    main_paths::Dict{String, String}
    model_dict::Union{ModelInputDictType, Nothing}
    materials_dict::Union{Dict{Any, Any}, Nothing}
    materials_library_dict::MaterialsDictType
    results_dict::Dict{String, Any}
    test_info_dict::Dict{String, Any}
    itime_step::Int64

    """ Benchmarks constructor

    # Inputs
    - `run_model::Bool`
        - Set to True to run model. Set to False to skip running model for
        cases where only post processing is desired.

    - `run_post_processing::Bool`
        - Set to True to run post processing.

    - `output_dir_root::Union{String, Nothing}`
        - Output directory where benchmark results will be sent.
    """
    function Benchmarks(;
        run_model::Bool = true,
        run_post_processing::Bool = true,
        output_dir_root::Union{String, Nothing} = nothing
    ) 
        # Use this if you want to run the benchmarks from the benchmarks directory
        # For this case place the models directory in the benchmarks directory
        project_path = dirname(@__DIR__)
        #project_path = pwd()
        println("Benchmarks project_path:", project_path)
        
        if isnothing(output_dir_root)
            output_dir_root = joinpath(
                project_path,
                "earthbox_benchmark_results_" * get_date_string()
            )
            println(">> output_dir_root:", output_dir_root)
        end

        check_output_dir_root(output_dir_root)

        lib = MaterialLibrary()
        materials_library_file = lib.benchmarks.path

        main_paths = make_main_paths_dict(
            project_path, output_dir_root, materials_library_file)

        new(
            run_model,
            run_post_processing,
            "",  # test_dir
            project_path,
            output_dir_root,
            materials_library_file,
            main_paths,
            nothing,  # model_dict
            nothing,  # materials_dict
            get_materials_library_dict(main_paths["materials_library_file"]),
            Dict{String, Any}(),
            Dict{String, Any}(),
            0
        )
    end
end

function get_date_string()
    return Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
end

function check_output_dir_root(output_dir_root::String)
    if !isdir(output_dir_root)
        mkpath(output_dir_root)
    end
end

function make_main_paths_dict(
    project_path::String,
    output_dir_root::String,
    materials_library_file::String
)::Dict{String, String}
    main_paths = Dict{String, String}()
    main_paths["materials_library_dir"] = dirname(materials_library_file)
    main_paths["materials_library_file"] = materials_library_file
    main_paths["project_path"] = project_path
    main_paths["models_path"] = joinpath(main_paths["project_path"], "models")
    main_paths["output_dir_root"] = output_dir_root
    main_paths["result_file_path"] = joinpath(
        output_dir_root, "test_results.yml"
    )
    return main_paths
end

function get_materials_library_dict(
    materials_library_file::String
)::MaterialsDictType
    (
        materials_library_dict, _
    ) = MaterialReader.read_materials_library(materials_library_file)
    return materials_library_dict
end

end # module

