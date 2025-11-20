module ModelParametersManager

Base.@kwdef mutable struct BlockGeometry
    xmin::Float64 = 0.4 # Normalized minimum x-coordinate of the block
    xmax::Float64 = 0.6 # Normalized maximum x-coordinate of the block
    ymin::Float64 = 0.4 # Normalized minimum y-coordinate of the block
    ymax::Float64 = 0.6 # Normalized maximum y-coordinate of the block
    zmin::Float64 = 0.4 # Normalized minimum z-coordinate of the block
    zmax::Float64 = 0.6 # Normalized maximum z-coordinate of the block
end

Base.@kwdef mutable struct BlockGeometry2d
    xmin::Float64 = 0.4 # Normalized minimum x-coordinate of the block
    xmax::Float64 = 0.6 # Normalized maximum x-coordinate of the block
    ymin::Float64 = 0.4 # Normalized minimum y-coordinate of the block
    ymax::Float64 = 0.6 # Normalized maximum y-coordinate of the block
end

function in_block(block_geometry::BlockGeometry, dxx::Float64, dyy::Float64, dzz::Float64)::Bool
    xmin = block_geometry.xmin
    xmax = block_geometry.xmax
    ymin = block_geometry.ymin
    ymax = block_geometry.ymax
    zmin = block_geometry.zmin
    zmax = block_geometry.zmax
    return xmin <= dxx <= xmax && ymin <= dyy <= ymax && zmin <= dzz <= zmax
end

function in_block(block_geometry::BlockGeometry2d, dxx::Float64, dyy::Float64)::Bool
    xmin = block_geometry.xmin
    xmax = block_geometry.xmax
    ymin = block_geometry.ymin
    ymax = block_geometry.ymax
    return xmin <= dxx <= xmax && ymin <= dyy <= ymax
end

Base.@kwdef mutable struct ModelParameters
    density_block::Float64 = 3100.0 # Block density, kg/m^3
    density_medium::Float64 = 3000.0 # Medium density, kg/m^3
    viscosity_block::Float64 = 1e+23 # Block viscosity, Pa*s
    viscosity_medium::Float64 = 1e+20 # Medium viscosity, Pa*s
    block_geometry::Union{BlockGeometry, BlockGeometry2d} = BlockGeometry()
end

function create_model_parameters(;model_type::Symbol=:ThreeDimensional)::ModelParameters
    model_parameters = ModelParameters()
    if model_type == :TwoDimensional
        model_parameters.block_geometry = BlockGeometry2d()
    end
    return model_parameters
end

end # module