module TtypeCalculator

import EarthBox.PrintFuncs: print_info
import EarthBox.DictUtils: check_for_missing_key_names
import EarthBox.EarthBoxDtypes: ParametersDictType

const REQUIRED_REFINEMENT_KEYS = [
    "xo_highres", "xf_highres", "dx_highres", "dx_lowres",
    "yf_highres", "dy_highres", "dy_lowres"
    ]

function get_minimum_ttype_refinement_parameters(
    type_refinement_parameters::Dict{String, Float64}
)
    xo_highres = get(type_refinement_parameters, "xo_highres", nothing)
    xf_highres = get(type_refinement_parameters, "xf_highres", nothing)
    dx_highres = get(type_refinement_parameters, "dx_highres", nothing)
    dx_lowres  = get(type_refinement_parameters, "dx_lowres", nothing)
    yf_highres = get(type_refinement_parameters, "yf_highres", nothing)
    dy_highres = get(type_refinement_parameters, "dy_highres", nothing)
    dy_lowres  = get(type_refinement_parameters, "dy_lowres", nothing)
    return xo_highres, xf_highres, dx_highres, dx_lowres, yf_highres, dy_highres, dy_lowres
end

""" Calculate the grid and marker parameters for type refinement.

(This is not currently used)

This function updates the number of basic nodes and markers per cell in the x 
and y directions for T-type grid refinement. The input dictionary 
`ttype_refinement_parameters` is checked to ensure that all required keys are present.

The number of markers in the x- and y-directions is updated only if the user provides
average marker spacings,`dx_marker` and `dy_marker` respectively. If the user does 
not provide these parameters, the number of markers per cell is not updated and
the values passed to the function are returned.

# Arguments
- `xsize::Float64`: Size of the model in the x direction in meters
- `ysize::Float64`: Size of the model in the y direction in meters
- `ttype_refinement_parameters::Dict{String, Float64}`: Dictionary of type refinement parameters
- `nmarkers_cell_x::Union{Float64, Nothing}`: Number of markers per cell in the x direction
- `nmarkers_cell_y::Union{Float64, Nothing}`: Number of markers per cell in the y direction
- `dx_marker::Union{Float64, Nothing}`: Spacing of markers in the x direction in meters
- `dy_marker::Union{Float64, Nothing}`: Spacing of markers in the y direction in meters

# Returns
- `xnum::Int64`: Number of basic nodes in the x direction
- `ynum::Int64`: Number of basic nodes in the y direction
- `nmarkers_cell_x::Union{Float64, Nothing}`: Number of markers per cell in the x direction
- `nmarkers_cell_y::Union{Float64, Nothing}`: Number of markers per cell in the y direction
"""
function calculate_grid_and_marker_parameters_for_type_refinement(
    xsize::Float64,
    ysize::Float64,
    ttype_refinement_parameters::Dict{String, Float64},
    nmarkers_cell_x::Union{Float64, Nothing}=nothing,
    nmarkers_cell_y::Union{Float64, Nothing}=nothing,
    dx_marker::Union{Float64, Nothing}=nothing,
    dy_marker::Union{Float64, Nothing}=nothing
)::Tuple{Int64, Int64, Union{Float64, Nothing}, Union{Float64, Nothing}}
    check_for_missing_key_names(REQUIRED_REFINEMENT_KEYS, ttype_refinement_parameters)
    xo_highres = ttype_refinement_parameters["xo_highres"]
    xf_highres = ttype_refinement_parameters["xf_highres"]
    dx_highres = ttype_refinement_parameters["dx_highres"]
    dx_lowres = ttype_refinement_parameters["dx_lowres"]
    yf_highres = ttype_refinement_parameters["yf_highres"]
    dy_highres = ttype_refinement_parameters["dy_highres"]
    dy_lowres = ttype_refinement_parameters["dy_lowres"]
    # Number of basic nodes in x and y directions
    xnum = calculate_number_of_basic_nodes_x(xsize, xo_highres, xf_highres, dx_highres, dx_lowres)
    ynum = calculate_number_of_basic_nodes_y(ysize, yf_highres, dy_highres, dy_lowres)
    # Average number of markers per cell in x and y directions
    if dx_marker !== nothing
        nmarkers_cell_x = calculate_number_of_markers_per_cell_x(
            xsize, xo_highres, xf_highres, dx_highres, dx_lowres, dx_marker)
    end
    if dy_marker !== nothing
        nmarkers_cell_y = calculate_number_of_markers_per_cell_y(
            ysize, yf_highres, dy_highres, dy_lowres, dy_marker)
    end
    return xnum, ynum, nmarkers_cell_x, nmarkers_cell_y
end

function ttype_parameters_provided_by_user(
    input_parameters::ParametersDictType
)::Bool
    xsize = get(input_parameters, "xsize", nothing)
    ysize = get(input_parameters, "ysize", nothing)
    xo_highres = get(input_parameters, "xo_highres", nothing)
    xf_highres = get(input_parameters, "xf_highres", nothing)
    dx_highres = get(input_parameters, "dx_highres", nothing)
    dx_lowres = get(input_parameters, "dx_lowres", nothing)
    yf_highres = get(input_parameters, "yf_highres", nothing)
    dy_highres = get(input_parameters, "dy_highres", nothing)
    dy_lowres = get(input_parameters, "dy_lowres", nothing)
    tmp = [xsize, ysize, xo_highres, xf_highres, dx_highres, dx_lowres, yf_highres, dy_highres, dy_lowres]
    none_count = count(==(nothing), tmp)
    check = false
    if none_count == 0
        check = true
    end
    return check
end

function calculate_basic_grid_parameters_for_type_refinement(
    input_parameters::ParametersDictType
)::Tuple{Int64, Int64}
    xsize = get(input_parameters, "xsize", nothing)
    ysize = get(input_parameters, "ysize", nothing)
    xo_highres = get(input_parameters, "xo_highres", nothing)
    xf_highres = get(input_parameters, "xf_highres", nothing)
    dx_highres = get(input_parameters, "dx_highres", nothing)
    dx_lowres = get(input_parameters, "dx_lowres", nothing)
    yf_highres = get(input_parameters, "yf_highres", nothing)
    dy_highres = get(input_parameters, "dy_highres", nothing)
    dy_lowres = get(input_parameters, "dy_lowres", nothing)
    tmp = [xsize, ysize, xo_highres, xf_highres, dx_highres, dx_lowres, yf_highres, dy_highres, dy_lowres]
    none_count = count(==(nothing), tmp)
    check = false
    if none_count > 0
        throw(ArgumentError("T-type parameters are incomplete. Fix your inputs."))
    else
        xnum = calculate_number_of_basic_nodes_x(xsize[1], xo_highres[1], xf_highres[1], dx_highres[1], dx_lowres[1])
        ynum = calculate_number_of_basic_nodes_y(ysize[1], yf_highres[1], dy_highres[1], dy_lowres[1])
        return xnum, ynum
    end
end

function calculate_number_of_basic_nodes_x(
    xsize::Float64,
    xo_highres::Float64,
    xf_highres::Float64,
    dx_highres::Float64,
    dx_lowres::Float64
)::Int64
    width_highres_x = calculate_width_high_res_x(xo_highres, xf_highres, xsize)
    width_lowres_x = calculate_width_low_res_x(xsize, xo_highres, xf_highres)
    number_of_basic_nodes_hrx = Int64(floor(width_highres_x/dx_highres))
    number_of_basic_nodes_lrx = Int64(floor(width_lowres_x/dx_lowres))
    result = number_of_basic_nodes_lrx*2 + number_of_basic_nodes_hrx + 1
    if result < 5
        throw(ArgumentError("The number of x basic nodes is less than 5. Decrease dx_lowres."))
    end
    print_info("Number of basic nodes x: $(result)", level=2)
    return result
end

function calculate_width_high_res_x(
    xo_highres::Float64, 
    xf_highres::Float64, 
    xsize::Float64
)::Float64
    result = xf_highres - xo_highres
    if result > xsize
        throw(ArgumentError(
            "The width of the high resolution region is greater than the model "
            *"width. Adjust xo_highres and xf_highres.")
            )
    end
    return result
end

function calculate_width_low_res_x(
    xsize::Float64,
    xo_highres::Float64,
    xf_highres::Float64,
)::Float64
    result = (xsize - calculate_width_high_res_x(xo_highres, xf_highres, xsize))/2.0
    if result < 0.0
        throw(ArgumentError(
            "The width of the low resolution region is less "
            *"than zero. Adjust xo_highres and xf_highres."
            )
        )
    end
    return result
end

function calculate_number_of_basic_nodes_y(
    ysize::Float64,
    yf_highres::Float64,
    dy_highres::Float64,
    dy_lowres::Float64
)::Int64
    number_of_basic_nodes_hr_y = Int64(floor(yf_highres/dy_highres))
    number_of_basic_nodes_lr_y = Int64(floor((ysize - yf_highres)/dy_lowres))
    result = number_of_basic_nodes_hr_y + number_of_basic_nodes_lr_y + 1
    if result < 5
        throw(ArgumentError("The number of y basic nodes is less than 5. Decrease dy_lowres."))
    end
    print_info("Number of basic nodes y: $(result)", level=2)
    return result
end

function calculate_number_of_markers_per_cell_x(
    xsize::Float64,
    xo_highres::Float64,
    xf_highres::Float64,
    dx_highres::Float64,
    dx_lowres::Float64,
    dx_marker::Float64
)::Float64
    xnum = calculate_number_of_basic_nodes_x(xo_highres, xf_highres, xsize, dx_highres, dx_lowres)
    avg_grid_spacing_x = xsize/(xnum - 1)
    result = Int64(floor(avg_grid_spacing_x/dx_marker))
    result = convert(Float64, result)
    if result < 2
        throw(ArgumentError("The number of x markers per cell is less than 2. Decrease dx_marker."))
    end
    print_info("Number of markers per cell x: $(result) with spacing $(dx_marker) m", level=1)
    return result
end

function calculate_number_of_markers_per_cell_y(
    ysize::Float64,
    yf_highres::Float64,
    dy_highres::Float64,
    dy_lowres::Float64,
    dy_marker::Float64
)::Float64
    ynum = calculate_number_of_basic_nodes_y(ysize, yf_highres, dy_highres, dy_lowres)
    avg_grid_spacing_y = ysize/(ynum - 1)
    result = Int64(floor(avg_grid_spacing_y/dy_marker))
    result = convert(Float64, result)
    if result < 2
        throw(ArgumentError("The number of y markers per cell is less than 2. Decrease dy_marker."))
    end
    print_info("Number of markers per cell y: $(result) with spacing $(dy_marker) m", level=1)
    return result
end

end # module

