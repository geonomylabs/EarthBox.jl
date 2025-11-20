module GridGeometryGroup

import EarthBox.Parameters: ParameterFloat, ParameterInt
import EarthBox.Parameters: print_parameter
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

mutable struct GridGeometry <: AbstractParameterGroup
    ynum::ParameterInt
    xnum::ParameterInt
    znum::ParameterInt
    ysize::ParameterFloat
    xsize::ParameterFloat
    zsize::ParameterFloat
    ystpavg::ParameterFloat
    xstpavg::ParameterFloat
    zstpavg::ParameterFloat
    ymin::ParameterFloat
    ymax::ParameterFloat
    xmin::ParameterFloat
    xmax::ParameterFloat
    zmin::ParameterFloat
    zmax::ParameterFloat
    xsize_start::ParameterFloat
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}

    function GridGeometry(
        ynum::ParameterInt,
        xnum::ParameterInt,
        znum::ParameterInt,
        ysize::ParameterFloat,
        xsize::ParameterFloat,
        zsize::ParameterFloat,
        ystpavg::ParameterFloat,
        xstpavg::ParameterFloat,
        zstpavg::ParameterFloat,
        ymin::ParameterFloat,
        ymax::ParameterFloat,
        xmin::ParameterFloat,
        xmax::ParameterFloat,
        zmin::ParameterFloat,
        zmax::ParameterFloat,
        xsize_start::ParameterFloat
    )
        obj_list = Union{ParameterFloat, ParameterInt}[]
        new(
            ynum, xnum, znum, ysize, xsize, zsize, ystpavg, xstpavg, zstpavg,
            ymin, ymax, xmin, xmax, zmin, zmax, xsize_start,
            obj_list
        )
    end
end

function GridGeometry(
    ynum::Int64,
    xnum::Int64,
    znum::Int64,
    ysize::Float64,
    xsize::Float64,
    zsize::Float64
)::GridGeometry
    data = GridGeometry(
        ParameterInt(ynum, "ynum", "None", "Basic grid resolution in y-direction."),
        ParameterInt(xnum, "xnum", "None", "Basic grid resolution in x-direction."),
        ParameterInt(znum, "znum", "None", "Basic grid resolution in z-direction."),
        ParameterFloat(ysize, "ysize", "m", "Height of model."),
        ParameterFloat(xsize, "xsize", "m", "Width of model."),
        ParameterFloat(zsize, "zsize", "m", "Depth of model."),
        ParameterFloat(0.0, "ystpavg", "m", "Average spacing of basic grid in y-direction."),
        ParameterFloat(0.0, "xstpavg", "m", "Average spacing of basic grid in x-direction."),
        ParameterFloat(0.0, "zstpavg", "m", "Average spacing of basic grid in z-direction."),
        ParameterFloat(0.0, "ymin", "m", "Minimum y-location"),
        ParameterFloat(ysize, "ymax", "m", "Maximum y-location"),
        ParameterFloat(0.0, "xmin", "m", "Minimum x-location"),
        ParameterFloat(xsize, "xmax", "m", "Maximum x-location"),
        ParameterFloat(0.0, "zmin", "m", "Minimum z-location"),
        ParameterFloat(zsize, "zmax", "m", "Maximum z-location"),
        ParameterFloat(xsize, "xsize_start", "m", "Starting value for xsize.")
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module 