"""
    ModelDataContainer

Container module with the constructor for the main model data structure and
the various model data structures contained within it.

"""
module ModelDataContainer

include("utils/ObjDictUtils.jl")
include("utils/DocTools.jl")
include("grids/Grids2dContainer.jl")
include("grids3d/Grids3dContainer.jl")
include("timestep/TimestepContainer.jl")
include("geometry/GeometryContainer.jl")
include("bcs/BoundaryConditionsContainer.jl")
include("markers/MarkerContainer.jl")
include("interpolation/InterpolationContainer.jl")
include("melting/MeltingContainer.jl")
include("topography/TopographyContainer.jl")
include("heat_equation/HeatEquationContainer.jl")
include("stokes_continuity/StokesContinuityContainer.jl")
include("multigrid/MultigridContainer.jl")
include("conversion/ConversionContainer.jl")
include("gravity/GravityContainer.jl")
include("benchmarks/BenchmarkContainer.jl")
include("carbonate/CarbonateContainer.jl")
include("materials/MaterialContainer.jl")
include("utils/OutputStandard.jl")

import EarthBox.Parameters: ParameterFloat
import EarthBox.Parameters: ParameterInt
import EarthBox.Parameters: ParameterStr
import EarthBox.Parameters: set_parameter_float!
import EarthBox.Parameters: set_parameter_int!
import EarthBox.Parameters: set_parameter_str!
import EarthBox.Parameters: convert_to_float64
import EarthBox.Parameters: convert_to_int64
import EarthBox.UnitConversion: UnitConversionData
import EarthBox.UnitConversion: get_standard_unit
import EarthBox.UnitConversion: get_conversion_func
import EarthBox.InputTools.Reader: get_parameters_input_dict
import EarthBox.ParameterRegistry: check_against_master
import EarthBox.EarthBoxDtypes: ParametersDictType
import EarthBox.EarthBoxDtypes: ObjDictType
import EarthBox.EarthBoxDtypes: GlobalIterDictType
import EarthBox: TtypeCalculator
import EarthBox.PrintFuncs: print_info
import .ObjDictUtils: add_collections_to_dict!
import .Grids2dContainer: Grids
import .Grids3dContainer: Grids3d
import .TimestepContainer: Timestep
import .GeometryContainer: Geometry
import .BoundaryConditionsContainer: BoundaryConditions
import .MarkerContainer: Markers
import .InterpolationContainer: Interpolation
import .MeltingContainer: Melting
import .TopographyContainer: Topography
import .HeatEquationContainer: HeatEquation
import .StokesContinuityContainer: StokesContinuity
import .GravityContainer: Gravity
import .BenchmarkContainer: Benchmarks
import .ConversionContainer: Conversion
import .CarbonateContainer: Carbonate
import .MaterialContainer: Materials

export ModelData, Grids, Grids3d, Timestep, Geometry, BoundaryConditions,
    Materials, Markers, Interpolation, Melting, Topography, HeatEquation, 
    StokesContinuity, Gravity, Benchmarks, Conversion, Carbonate, RGAS

"""
    ModelData

Struct containing parameter and array collections for the EarthBox model data and other fields 
related to the EarthBox model. Each collection container contains parameter and array collections, 
`parameters` and `arrays` respectively, for the various data categories. Parameters and array 
collections are in turn composed of groups of related parameters and arrays. The general structure 
of the model data collections container is as follows:

```
ModelData
├─ collection_container_1
│  ├─ parameters
│  │  ├─ parameter_group_A
│  │  │  └─ parameter_A1
│  │  │  └─ parameter_A2
│  │  └─ parameter_group_B
│  │     ├─ parameter_B1
│  │     ├─ parameter_B2
│  └─ arrays
│     ├─ array_group_C
│     │  └─ array_C1
│     │  └─ array_C2
│     └─ array_group_D
│        └─ array_D1
│        └─ array_D2
├─ collection_container_2
│  ├─ parameters
│  │  ├─ parameter_group_E
│  │  │  └─ parameter_E1
│  │  │  └─ parameter_E2
│  │  └─ parameter_group_F
│  │     ├─ parameter_F1
│  │     ├─ parameter_F2
│  └─ arrays
│     ├─ array_group_G
│     │  └─ array_G1
│     │  └─ array_G2
│     └─ array_group_H
│        └─ array_H1
│        └─ array_H2
```

For a ModelData object `model`, the `parameter_A1` and `array_C1` can be accessed as follows: 
- `model.collection_container_1.parameters.parameter_group_A.parameter_A1`
- `model.collection_container_1.arrays.array_group_C.array_C1`

Parameter groups contain parameter objects with attributes like `value`, `unit`, and `description`.
Array groups contain array objects with attributes like `array`, `units`, and `description`.

# Fields
- `grids::Grids` or `Grids3d`: Parameter and array collections for basic and staggered grids. See [`Grids`](@ref) for details.
- `timestep::Timestep`: Parameter collections for model time stepping. See [`Timestep`](@ref) for details.
- `geometry::Geometry`: Parameter collections for material geometry. See [`Geometry`](@ref) for details.
- `bcs::BoundaryConditions`: Parameter and array collections for boundary conditions. See [`BoundaryConditions`](@ref) for details.
- `materials::Materials`: Parameter, array and dictionary collections for material properties. See [`Materials`](@ref) for details.
- `markers::Markers`: Parameter and array collections for markers. See [`Markers`](@ref) for details.
- `interpolation::Interpolation`: Parameter and array collections for marker-grid interpolation. See [`Interpolation`](@ref) for details.
- `melting::Melting`: Parameter and array collections for melting and melt extraction. See [`Melting`](@ref) for details.
- `topography::Topography`: Parameter and array collections for topography evolution. See [`Topography`](@ref) for details.
- `heat_equation::HeatEquation`: Parameter and array collections for heat equation. See [`HeatEquation`](@ref) for details.
- `stokes_continuity::StokesContinuity`: Parameter and array collections for Stokes-continuity solver. See [`StokesContinuity`](@ref) for details.
- `gravity`: `Gravity` parameters including gravity components and control flags. See [`Gravity`](@ref) for details.
- `benchmarks::Benchmarks`: Parameter collection for benchmark parameters. See [`Benchmarks`](@ref) for details.
- `conversion::Conversion`: Conversion parameters. See [`Conversion`](@ref) for details.
- `carbonate::Carbonate`: Parameter and array collections for carbonate deposition configuration. See [`Carbonate`](@ref) for details.
- `RGAS::ParameterFloat`: Gas constant. See [`ParameterFloat`](@ref) for details on type.
- `unit_conversion_data::UnitConversionData`: Data structure with unit conversion functions and registry.
- `model_input_file::Union{String, Nothing}`: Model input file.
- `input_parameters::ParametersDictType`: Input parameters dictionary.
- `obj_dict::ObjDictType`: Dictionary of EarthBox parameter and array objects where keys are names and values are parameter and array objects.
- `global_iter_dict::GlobalIterDictType`: Dictionary of global iteration information.

"""
mutable struct ModelData
    grids::Union{Grids, Grids3d}
    timestep::Timestep
    geometry::Geometry
    bcs::BoundaryConditions
    materials::Materials
    markers::Markers
    interpolation::Interpolation
    melting::Melting
    topography::Topography
    heat_equation::HeatEquation
    stokes_continuity::StokesContinuity
    gravity::Gravity
    benchmarks::Benchmarks
    conversion::Conversion
    carbonate::Carbonate
    RGAS::ParameterFloat
    unit_conversion_data::UnitConversionData
    model_input_file::Union{String, Nothing}
    input_parameters::ParametersDictType
    obj_dict::ObjDictType
    global_iter_dict::GlobalIterDictType

    function ModelData(
        model_input_file::Union{String, Nothing},
        initialization_params::Union{ParametersDictType, Nothing}
    )::ModelData
        unit_conversion_data = UnitConversionData()
        input_parameters = get_input_parameters_dict(model_input_file, initialization_params)
        ynum, xnum, znum = get_nums(input_parameters)
        update_input_parameters_with_nums!(input_parameters, ynum, xnum)
        xsize, ysize, zsize = get_xsize_and_ysize_in_model_units(input_parameters)
        nmarkers_cell_y, nmarkers_cell_x, nmarkers_cell_z = get_nmarkers(input_parameters)
        dx_marker, dy_marker = get_marker_spacing(input_parameters)
        
        marker_parameters = calculate_marker_parameters_tuple(
            ynum, xnum, znum, ysize, xsize, zsize, 
            nmarkers_cell_y, nmarkers_cell_x, nmarkers_cell_z,
            dx_marker, dy_marker
            )

        _topography = Topography()
        toponum = _topography.parameters.topo_grid.toponum.value

        model = new(
            get_grids(ynum, xnum, znum, ysize, xsize, zsize),          # grids
            Timestep(),                                                # timestep
            Geometry(),                                                # geometry
            BoundaryConditions(ynum, xnum),                            # bcs
            Materials(),                                               # materials
            Markers(marker_parameters),                                # markers
            Interpolation(ynum, xnum, marker_parameters.marknum),      # interpolation
            Melting(),                                                 # melting
            _topography,                                               # topography
            HeatEquation(ynum, xnum),                                  # heat_equation
            StokesContinuity(ynum, xnum),                              # stokes_continuity
            Gravity(),                                                 # gravity
            Benchmarks(),                                              # benchmarks
            Conversion(),                                              # conversion
            Carbonate(toponum),                                        # carbonate
            ParameterFloat(8.314, "RGAS", "J/mol/K", "Gas constant."), # RGAS
            unit_conversion_data,                                      # unit_conversion_data
            model_input_file,                                          # model_input_file
            input_parameters,                                          # input_parameters
            ObjDictType(),                                             # obj_dict
            GlobalIterDictType()                                       # global_iter_dict
        )
        
        add_objects_to_dict!(model)
        update_parameter_units!(model)
        load_input_parameters!(model)
        
        return model
    end
end

function get_nums(input_parameters::ParametersDictType)::Tuple{Int64, Int64, Int64}
    if TtypeCalculator.ttype_parameters_provided_by_user(input_parameters)
        print_info("Number of basic nodes in x and y directions are calculated from T-type refinement parameters", level=1)
        xnum, ynum = TtypeCalculator.calculate_basic_grid_parameters_for_type_refinement(input_parameters)
    else
        print_info("Number of basic nodes in x and y directions are provided as input parameters", level=1)
        ynum = input_parameters["ynum"][1]
        xnum = input_parameters["xnum"][1]
    end
    znum = get(input_parameters, "znum", [1])[1]
    @assert xnum > 3 "Number of basic nodes in x direction must be greater than 1."
    @assert ynum > 3 "Number of basic nodes in y direction must be greater than 1."
    return (ynum, xnum, znum)
end

function update_input_parameters_with_nums!(
    input_parameters::ParametersDictType,
    ynum::Int64,
    xnum::Int64
)::Nothing
    @assert xnum > 3 "Number of basic nodes in x direction must be greater than 1."
    @assert ynum > 3 "Number of basic nodes in y direction must be greater than 1."
    input_parameters["xnum"] = [xnum, "None"]
    input_parameters["ynum"] = [ynum, "None"]
    return nothing
end

function get_nmarkers(input_parameters::ParametersDictType)::Tuple{Float64, Float64, Float64}
    if !marker_xspacing_provided_by_user(input_parameters)
        print_info("Number of markers per cell in x direction are provided as input parameters", level=1)
        nmarkers_cell_x = input_parameters["nmarkers_cell_x"][1]
        @assert nmarkers_cell_x > 0.0 "Number of markers per cell in x direction must be greater than 0."
    else
        nmarkers_cell_x = calculate_number_of_markers_per_cell(
            input_parameters["xsize"][1],
            input_parameters["xnum"][1],
            input_parameters["dx_marker"][1]
        )
        print_info("Number of calculated markers per cell in x direction: $(nmarkers_cell_x)", level=1)
        @assert nmarkers_cell_x > 0.0 "Number of markers per cell in x direction must be greater than 0."
    end
    if !marker_yspacing_provided_by_user(input_parameters)
        print_info("Number of markers per cell in y direction are provided as input parameters", level=1)
        nmarkers_cell_y = input_parameters["nmarkers_cell_y"][1]
        @assert nmarkers_cell_y > 0.0 "Number of markers per cell in y direction must be greater than 0."
    else
        nmarkers_cell_y = calculate_number_of_markers_per_cell(
            input_parameters["ysize"][1],
            input_parameters["ynum"][1],
            input_parameters["dy_marker"][1]
        )
        print_info("Number of calculated markers per cell in y-direction : $(nmarkers_cell_y)", level=1)
        @assert nmarkers_cell_y > 0.0 "Number of markers per cell in y direction must be greater than 0."
    end
    nmarkers_cell_z = get(input_parameters, "nmarkers_cell_z", [1.0])[1]
    return (nmarkers_cell_y, nmarkers_cell_x, nmarkers_cell_z)
end

function marker_xspacing_provided_by_user(input_parameters::ParametersDictType)::Bool
    dx_marker_info = get(input_parameters, "dx_marker", nothing)
    return dx_marker_info !== nothing
end

function marker_yspacing_provided_by_user(input_parameters::ParametersDictType)::Bool
    dy_marker_info = get(input_parameters, "dy_marker", nothing)
    return dy_marker_info !== nothing
end

function get_marker_spacing(
    input_parameters::ParametersDictType
)::Tuple{Union{Float64, Nothing}, Union{Float64, Nothing}}
    dx_marker_info = get(input_parameters, "dx_marker", nothing)
    if dx_marker_info !== nothing
        dx_marker = dx_marker_info[1]
    else
        dx_marker = nothing
    end
    dy_marker_info = get(input_parameters, "dy_marker", nothing)
    if dy_marker_info !== nothing
        dy_marker = dy_marker_info[1]
    else
        dy_marker = nothing
    end
    return dx_marker, dy_marker
end

function calculate_number_of_markers_per_cell(
    length::Float64,
    number_of_nodes::Int64,
    marker_spacing::Float64
)::Float64
    avg_grid_spacing = length/(number_of_nodes - 1)
    nmarkers_per_cell = Int64(floor(avg_grid_spacing/marker_spacing))
    nmarkers_per_cell = convert(Float64, nmarkers_per_cell)
    return nmarkers_per_cell
end

function get_grids(
    ynum::Int64, xnum::Int64, znum::Int64, 
    ysize::Float64, xsize::Float64, zsize::Float64
)::Union{Grids, Grids3d}
    if znum > 1 && zsize > 0.0
        return Grids3d(ynum=ynum, xnum=xnum, znum=znum, ysize=ysize, xsize=xsize, zsize=zsize)
    else
        return Grids(ynum, xnum, ysize, xsize)
    end
end

function calculate_marker_parameters_tuple(
    ynum::Int,
    xnum::Int,
    znum::Int,
    ysize::Float64,
    xsize::Float64,
    zsize::Float64,
    nmarkers_cell_y::Float64,
    nmarkers_cell_x::Float64,
    nmarkers_cell_z::Float64,
    dx_marker::Union{Float64, Nothing},
    dy_marker::Union{Float64, Nothing}
)::NamedTuple
    mxnum, mynum, mznum = calculate_number_of_markers_for_each_dimension(
        xnum, ynum, znum, nmarkers_cell_x, nmarkers_cell_y, nmarkers_cell_z
    )
    marknum = calculate_total_number_of_markers(mxnum, mynum, mznum)
    dx_marker, dy_marker = calculate_marker_spacing(
        xsize, ysize, xnum, ynum, nmarkers_cell_x, 
        nmarkers_cell_y, dx_marker, dy_marker
    )
    @assert dx_marker > 0.0 "Marker spacing in x direction must be positive."
    @assert dy_marker > 0.0 "Marker spacing in y direction must be positive."
    return (
        nmarkers_cell_y = nmarkers_cell_y,
        nmarkers_cell_x = nmarkers_cell_x,
        nmarkers_cell_z = nmarkers_cell_z,
        mynum = mynum,
        mxnum = mxnum,
        mznum = mznum,
        marknum = marknum,
        mystep = ysize / mynum,
        mxstep = xsize / mxnum,
        mzstep = zsize / mznum,
        dx_marker = dx_marker,
        dy_marker = dy_marker
    )
end

function calculate_number_of_markers_for_each_dimension(
    xnum::Int64,
    ynum::Int64,
    znum::Int64,
    nmarkers_cell_x::Float64,
    nmarkers_cell_y::Float64,
    nmarkers_cell_z::Float64
)
    mxnum = Int((xnum - 1) * nmarkers_cell_x)
    mynum = Int((ynum - 1) * nmarkers_cell_y)
    mznum = Int((znum - 1) * nmarkers_cell_z)
    if mznum < 1
        mznum = 1
    end
    return mxnum, mynum, mznum
end

function calculate_total_number_of_markers(
    mxnum::Int64,
    mynum::Int64,
    mznum::Int64
)::Int64
    return mxnum * mynum * mznum
end

function calculate_marker_spacing(
    xsize::Float64,
    ysize::Float64,
    xnum::Int64,
    ynum::Int64,
    nmarkers_cell_x::Float64,
    nmarkers_cell_y::Float64,
    dx_marker::Union{Float64, Nothing},
    dy_marker::Union{Float64, Nothing}
)::Tuple{Union{Float64, Nothing}, Union{Float64, Nothing}}
    # if marker spacing was not provided as an input, calculate it
    if dx_marker === nothing
        dx_marker = xsize / (nmarkers_cell_x*(xnum-1))
    end
    if dy_marker === nothing
        dy_marker = ysize / (nmarkers_cell_y*(ynum-1))
    end
    return dx_marker, dy_marker
end

""" Get input parameters dictionary

The input_parameters dictionary contains either key input for model initialization or all input
parameters read from the model input file. 

The input_parameters dictionary will be set equal to the initialization_params dictionary if defined. 
Otherwise, the input_parameters dictionary will be set equal to the input parameters read from the 
model input file.  The keys of the initialize_params dictionary are dependent on the inputs 
provides by the user via the API. For example, if the user provides input for a T-type refined grid
the input_parameters dictionary will contain the keys "xsize", "ysize", "xo_highres", "xf_highres", 
"dx_highres", "yf_highres", and "dy_highres" whereas if the user provides input for a general basic 
grid the input_parameters dictionary will contain the keys "xsize", "ysize", "xnum", and "ynum".
EarthBox will calculate the remaining input parameters based on the user provided inputs.

The initialization_params dictionary is a collection parameters for initializing grid and marker
dimensions and resolution provided as input in the API. Additional input parameter must be defined
by the user during subsequent API calls.

If a model input file is provided, all input parameters will be read from the input file.

"""
function get_input_parameters_dict(
    model_input_file::Union{String, Nothing},
    initialization_params::Union{ParametersDictType, Nothing}
)::ParametersDictType
    input_parameters = nothing
    if initialization_params !== nothing
        print_info(
            "Key initialization input parameters are provided from the API. "
            *"Additional parameters must be entered via subsequent API calls."
            )
        input_parameters = initialization_params
    elseif model_input_file !== nothing
        print_info("Key initialization input parameters are read from the model input file.")
        input_parameters = get_parameters_input_dict(model_input_file)
    else
        throw(
            ArgumentError(
            "The initialization_params dictionary and model_input_file are set to nothing. Check your inputs!."
            )
            )
    end
    return input_parameters
end

function add_objects_to_dict!(model::ModelData)::Nothing
    add_collections_to_dict!(model.obj_dict, model.grids)
    add_collections_to_dict!(model.obj_dict, model.timestep)
    add_collections_to_dict!(model.obj_dict, model.geometry)
    add_collections_to_dict!(model.obj_dict, model.bcs)
    add_collections_to_dict!(model.obj_dict, model.materials)
    add_collections_to_dict!(model.obj_dict, model.markers)
    add_collections_to_dict!(model.obj_dict, model.interpolation)
    add_collections_to_dict!(model.obj_dict, model.melting)
    add_collections_to_dict!(model.obj_dict, model.topography)
    add_collections_to_dict!(model.obj_dict, model.heat_equation)
    add_collections_to_dict!(model.obj_dict, model.stokes_continuity)
    add_collections_to_dict!(model.obj_dict, model.gravity)
    add_collections_to_dict!(model.obj_dict, model.benchmarks)
    add_collections_to_dict!(model.obj_dict, model.conversion)
    add_collections_to_dict!(model.obj_dict, model.carbonate)
    return nothing
end

function update_parameter_units!(model::ModelData)::Nothing
    for (parameter_name, parameter_list) in model.input_parameters
        param_obj = get_object(model, parameter_name)
        unit_start_init = parameter_list[2]
        unit_end_init = param_obj.units
        conversion_func, unit_start, unit_end = get_conversion_func(
            unit_start_init, unit_end_init, model.unit_conversion_data, parameter_name)
        value_start = parameter_list[1]
        value_end = conversion_func(value_start)
        update_parameter!(model, parameter_name, unit_end, value_end)
        if unit_start != unit_end
            print_info("$parameter_name : $value_start $unit_start -> $value_end $unit_end", level=2)
        end
    end
    return nothing
end

function update_parameter!(
    model::ModelData,
    name::String,
    units::String,
    value::Union{Float64, Int64, String}
)::Nothing
    model.input_parameters[name][1] = value
    model.input_parameters[name][2] = units
    return nothing
end

function load_input_parameters!(model::ModelData)::Nothing
    print_info("Loading input parameters", level=1)
    for (name, info) in model.input_parameters
        value = info[1]
        units = info[2]
        units = get_standard_unit(units, model.unit_conversion_data.unit_registry)
        obj = get_object(model, name)
        obj_units = obj.units
        obj_units = get_standard_unit(obj_units, model.unit_conversion_data.unit_registry)
        
        if units != obj_units
            println("name: $name")
            println("units: $units")
            println("obj_units: $obj_units")
            throw(ArgumentError(
                "Input parameter unit $units is not equal to the object unit "
                *"$obj_units for parameter $name."
            ))
        end
        
        set_parameter!(model, name, value)
        print_info("$name with value $value and units $units", level=2)
    end
    return nothing
end

function get_object(model::ModelData, obj_name::String)::Any
    if !haskey(model.obj_dict, obj_name)
        throw(ArgumentError(
            "Parameter name $obj_name is not in obj_dict. "
            *"Check that the parameter name is correct and."
        ))
    end
    return model.obj_dict[obj_name]
end

function set_model_data!(
    model::ModelData,
    parameters::Union{Dict{String, Union{Float64, Int64}}, Nothing}=nothing,
    material_domain_ids::Union{Dict{String, Int64}, Nothing}=nothing
)::Nothing
    if parameters !== nothing
        set_parameters!(model, parameters)
    end
    if material_domain_ids !== nothing
        set_material_domain_ids!(model, material_domain_ids)
    end
    return nothing
end

function load_parameters!(
    model::ModelData,
    valid_names::Tuple{Vararg{Symbol}};
    kwargs...
)::Nothing
    for (key, value) in kwargs
        if (key in valid_names)
            set_parameter!(model, String(key), value)
        else
            print_valid_object_info(model, valid_names)
            error("Invalid input parameter name: $key")
        end
    end
    return nothing
end

function print_valid_object_info(
    model::ModelData, 
    valid_names::Tuple{Vararg{Symbol}}
)::Nothing
    println("Valid parameters are as follows:")
    for symbol_name in valid_names
        name = String(symbol_name)
        obj = get_object(model, name)
        println(obj.name, " (", obj.units, "): ", obj.description)
    end
    return nothing
end

function set_parameters!(
    model::ModelData,
    new_parameters::Union{Dict{String, <:Union{Float64, Int64, String}}, Nothing}=nothing
)::Nothing
    if new_parameters !== nothing
        for (param_name, value) in new_parameters
            set_parameter!(model, param_name, value)
        end
    end
    return nothing
end

function set_parameter!(
    model::ModelData,
    param_name::Union{Symbol, String},
    value::Union{Float64, Int64, String, Bool, Nothing}
)::Nothing
    param_name = String(param_name)
    param_obj = get_object(model, param_name)
    if param_obj isa ParameterFloat && value !== nothing
        set_parameter_float!(param_obj, convert_to_float64(value))
    elseif param_obj isa ParameterInt && value !== nothing
        set_parameter_int!(param_obj, convert_to_int64(value))
    elseif param_obj isa ParameterStr && value !== nothing
        set_parameter_str!(param_obj, string(value))
    end
    return nothing
end

function to_float64(value::Union{Float64, Int64, String, Bool})::Float64
    if typeof(value) == Float64
        return value
    elseif typeof(value) == Int64
        return convert(Float64, value)
    elseif typeof(value) == String
        return parse(Float64, value)
    elseif typeof(value) == Bool
        return value ? 1.0 : 0.0
    end
end

function set_material_domain_ids!(
    model::ModelData,
    new_matid_domains::Dict{String, Int64}
)::Nothing
    for (domain_name, matid) in new_matid_domains
        try
            model.materials.dicts.matid_domains[domain_name] = matid
        catch e
            throw(ArgumentError("Material domain name $domain_name is not valid."))
        end
    end
    return nothing
end

function check_for_name(iflag::Int64, target_name::String, obj_dict::ObjDictType)::Nothing
    keys_list = collect(keys(obj_dict))
    if target_name in keys_list
        println("found target name $target_name for flag $iflag")
    end
    return nothing
end

function get_xsize_and_ysize_in_model_units(
    parameters::ParametersDictType
)::Tuple{Float64, Float64, Float64}
    xsize = convert_dimension_to_model_units(parameters, "xsize")
    ysize = convert_dimension_to_model_units(parameters, "ysize")
    zsize = convert_dimension_to_model_units(parameters, "zsize")
    @assert xsize > 0.0 "Model size in x direction must be positive."
    @assert ysize > 0.0 "Model size in y direction must be positive."
    return (xsize, ysize, zsize)
end

function convert_dimension_to_model_units(
    parameters::ParametersDictType, 
    dimension_name::String
)::Float64
    if !haskey(parameters, dimension_name)
        return 0.0
    end
    dimension = parameters[dimension_name][1]
    unit_start = parameters[dimension_name][2]
    unit_end = "m"
    unit_conversion_data = UnitConversionData()
    conversion_func, unit_start, unit_end = get_conversion_func(
        unit_start, unit_end, unit_conversion_data, dimension_name)
    return conversion_func(dimension)
end

end # module 