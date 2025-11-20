module Convergence

using YAML
import EarthBox.EarthBoxDtypes: GlobalIterDictType

function save_iteration_info_to_file(
    output_dir::String,
    iter_dict::GlobalIterDictType
)::Nothing
    file_path = joinpath(output_dir, "convergence.yml")
    make_yaml_file(file_path, iter_dict)
    return nothing
end

function make_yaml_file(yaml_file_path::String, data_dict)
    open(yaml_file_path, "w") do ymlfile
        YAML.write(ymlfile, data_dict)
    end
end

end # module 