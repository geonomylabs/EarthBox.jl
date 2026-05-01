""" Make initialization parameter dictionary

The initialization_params dictionary provides the user a way to define key
parameters for initializing the model via the API including model dimensions,
nodal resolution and grid refinement parameters. If this information is not
provided (i.e. initialization_params = nothing), the model will be initialized
with parameters read from the model input file.

"""
function make_initialization_parameter_dict(
    xnum::Union{Int, Nothing},
    ynum::Union{Int, Nothing},
    znum::Union{Int, Nothing},
    dx_marker::Union{Float64, Nothing},
    dy_marker::Union{Float64, Nothing},
    dz_marker::Union{Float64, Nothing},
    nmarkers_cell_x::Union{Float64, Nothing},
    nmarkers_cell_y::Union{Float64, Nothing},
    nmarkers_cell_z::Union{Float64, Nothing},
    xsize::Union{Float64, Nothing},
    ysize::Union{Float64, Nothing},
    zsize::Union{Float64, Nothing},
    type_refinement_parameters::Union{Dict{String, Float64}, Nothing}
)::Union{ParametersDictType, Nothing}
    if type_refinement_parameters !== nothing
        initialization_params = define_init_parameters_for_ttype_grid(
            xsize, ysize, type_refinement_parameters
        )
    else
        initialization_params = define_init_parameters_for_general_basic_grid(
            xsize, ysize, zsize, xnum, ynum, znum
        )
    end
    add_marker_resolution_parameters!(
        initialization_params,
        dx_marker, dy_marker, dz_marker,
        nmarkers_cell_x, nmarkers_cell_y, nmarkers_cell_z
    )
    return initialization_params
end

function define_init_parameters_for_general_basic_grid(
    xsize::Union{Float64, Nothing},
    ysize::Union{Float64, Nothing},
    zsize::Union{Float64, Nothing},
    xnum::Union{Int, Nothing},
    ynum::Union{Int, Nothing},
    znum::Union{Int, Nothing}
)::Union{ParametersDictType, Nothing}
    # These minimum attributes are used to determine if the user has provided
    # the minimum number of input parameters to initialize the model. If the
    # has not provided the minimum number of input parameters nothing is returned.
    # This instructs the code to use parameters from a model input file as
    # opposed to using the user provided parameters via the API.
    minimum_attr_list = [xsize, ysize, xnum, ynum]
    none_count = count(==(nothing), minimum_attr_list)
    # If some but not all key initialization parameters are provided, throw an error as the user
    # has not provided sufficient inputs for initializing the model and ambiguity is introduced.
    if 0 < none_count < length(minimum_attr_list)
        throw(ArgumentError(
            "Sufficient API inputs were not provided for initializing the general basic grid:  "
            *"xsize = $xsize, ysize = $ysize, xnum = $xnum, ynum = $ynum. Either provide all "
            *"inputs via the API or define entries for these parameters in the model input file."
            )
            )
    end
    # If all key initialization parameters are provided (none_count = 0), return the
    # initialization_params dictionary
    if none_count == 0
        print_info("Sufficient API inputs were provided for initializing the general basic grid.")
        initialization_params =  ParametersDictType(
            "xnum" => [Int(xnum), "None"],
            "ynum" => [Int(ynum), "None"],
            "xsize" => [Float64(xsize), "m"],
            "ysize" => [Float64(ysize), "m"]
        )
        # Add 3d parameters if they are provided by the user.
        if znum !== nothing
            initialization_params["znum"] = [Int(znum), "None"]
        end
        if zsize !== nothing
            initialization_params["zsize"] = [Float64(zsize), "m"]
        end
        return initialization_params
    else
        return nothing
    end
end

function define_init_parameters_for_ttype_grid(
    xsize::Union{Float64, Nothing},
    ysize::Union{Float64, Nothing},
    type_refinement_parameters::Dict{String, Float64}
)::Union{ParametersDictType, Nothing}
    (
        xo_highres, xf_highres, dx_highres, dx_lowres,
        yf_highres, dy_highres, dy_lowres
    ) = TtypeCalculator.get_minimum_ttype_refinement_parameters(type_refinement_parameters)
    minimum_attr_list = [
        xsize, ysize,
        xo_highres, xf_highres, dx_highres, dx_lowres,
        yf_highres, dy_highres, dy_lowres
    ]
    # These minimum attributes are used to determine if the user has provided
    # the minimum number of input parameters to initialize the model. If the
    # has not provided the minimum number of input parameters nothing is returned.
    # This instructs the code to use parameters from a model input file as
    # opposed to using the user provided parameters via the API.
    none_count = count(==(nothing), minimum_attr_list)
    if 0 < none_count < length(minimum_attr_list)
        throw(ArgumentError(
            "Sufficient API inputs were not provided for initializing the T-type refined grid:  "
            *"xsize = $xsize, ysize = $ysize, xo_highres = $xo_highres, xf_highres = $xf_highres, "
            *"dx_highres = $dx_highres, dx_lowres = $dx_lowres, yf_highres = $yf_highres, "
            *"dy_highres = $dy_highres, dy_lowres = $dy_lowres. Either provide all "
            *"inputs via the API or define entries for these parameters in the model input file."
            )
            )
    end
    if none_count == 0
        print_info(
            "Parameters for initializing the T-type refined grid are provided as "
            * "input parameters from the API"
            )
        initialization_params =  ParametersDictType(
            "xsize" => [Float64(xsize), "m"],
            "ysize" => [Float64(ysize), "m"],
            "xo_highres" => [Float64(xo_highres), "m"],
            "xf_highres" => [Float64(xf_highres), "m"],
            "dx_highres" => [Float64(dx_highres), "m"],
            "dx_lowres" => [Float64(dx_lowres), "m"],
            "yf_highres" => [Float64(yf_highres), "m"],
            "dy_highres" => [Float64(dy_highres), "m"],
            "dy_lowres" => [Float64(dy_lowres), "m"]
        )
        return initialization_params
    else
        return nothing
    end
end

function add_marker_resolution_parameters!(
    initialization_params::Union{ParametersDictType, Nothing},
    dx_marker::Union{Float64, Nothing},
    dy_marker::Union{Float64, Nothing},
    dz_marker::Union{Float64, Nothing},
    nmarkers_cell_x::Union{Float64, Nothing},
    nmarkers_cell_y::Union{Float64, Nothing},
    nmarkers_cell_z::Union{Float64, Nothing}
)::Nothing

    if dx_marker !== nothing && initialization_params !== nothing
        print_info("Marker spacing in x direction is provided as input parameter from the API")
        initialization_params["dx_marker"] = [Float64(dx_marker), "m"]
    elseif nmarkers_cell_x !== nothing && initialization_params !== nothing
        print_info("Number of markers per cell in x direction is provided as input parameter from the API")
        initialization_params["nmarkers_cell_x"] = [Float64(nmarkers_cell_x), "None"]
    elseif nmarkers_cell_x == nothing && dx_marker == nothing && initialization_params !== nothing
        error(
            "You are trying to initialize the model through the API but have not provided either the "
            * "marker spacing in x direction (dx_marker) or the number of markers per cell in the "
            * "x-direction (nmarkers_cell_x). One of these parameters are required to initialize the "
            * "marker swarm."
            )
    end

    if dy_marker !== nothing && initialization_params !== nothing
        print_info("Marker spacing in y direction is provided as input parameter from the API")
        initialization_params["dy_marker"] = [Float64(dy_marker), "m"]
    elseif nmarkers_cell_y !== nothing && initialization_params !== nothing
        print_info("Number of markers per cell in y direction is provided as input parameter from the API")
        initialization_params["nmarkers_cell_y"] = [Float64(nmarkers_cell_y), "None"]
    elseif nmarkers_cell_y == nothing && dy_marker == nothing && initialization_params !== nothing
        error(
            "You are trying to initialize the model through the API but have not provided either the "
            * "marker spacing in y direction (dy_marker) or the number of markers per cell in the "
            * "y-direction (nmarkers_cell_y). One of these parameters are required to initialize the "
            * "marker swarm."
            )
    end
    # Add statement for 3D once implemented
    #if dz_marker !== nothing
    #    initialization_params["dz_marker"] = [Float64(dz_marker), "m"]
    #else
    #    initialization_params["nmarkers_cell_z"] = [Float64(nmarkers_cell_z), "None"]
    #end
    return nothing
end

function check_inputs(
    eb_paths::EarthBoxPaths.EarthBoxPathsState,
    initialization_params::Union{ParametersDictType, Nothing}
)::Nothing
    if !haskey(eb_paths.paths, "model_input_file") && initialization_params === nothing
        throw(ArgumentError(
            "Both `model_input_file` and initialization parameters are not defined. EarthBox "
            *"cannot be initialized. One of these must be defined."
            ))
    end
    if haskey(eb_paths.paths, "model_input_file") && initialization_params !== nothing
        print_warning(
            "Model input file is provided and key initialization parameters are provided as input "
            *"parameters from the API. The model input file will be ignored."
            )
    end
    return nothing
end
