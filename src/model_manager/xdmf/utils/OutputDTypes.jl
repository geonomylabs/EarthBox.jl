module OutputDTypes

""" Struct used to store jld2 information associated with mesh.
"""
struct Mesh2djld
    xnum::Int
    ynum::Int
    ysize::Float64
    xsize::Float64
    jld_mesh_filename::String
    jld_dataname_1dy::String
    jld_dataname_1dx::String
    jld_dataname_pr1dy::String
    jld_dataname_pr1dx::String
    gridy_km_array::Vector{Float64}
    gridx_km_array::Vector{Float64}
    gridy_pr_km_array::Vector{Float64}
    gridx_pr_km_array::Vector{Float64}
end

function Mesh2djld(mesh_data::Dict{String, Any})
    return Mesh2djld(
        mesh_data["xnum"],
        mesh_data["ynum"],
        mesh_data["ysize"],
        mesh_data["xsize"],
        mesh_data["jld_mesh_filename"],
        mesh_data["jld_dataname_1dy"],
        mesh_data["jld_dataname_1dx"],
        mesh_data["jld_dataname_pr1dy"],
        mesh_data["jld_dataname_pr1dx"],
        mesh_data["gridy_km_array"],
        mesh_data["gridx_km_array"],
        mesh_data["gridy_pr_km_array"],
        mesh_data["gridx_pr_km_array"]
    )
end

""" Struct used to define a scalar field within a jld2 file.
"""
struct ScalarField
    name::String
    jld_dataname::String
    number_type::String
    units::String
    header::String
    grid_type::String
    scalar_array::Union{Vector{Float64}, Vector{Int64}, Matrix{Float64}, Matrix{Int64}}
end

""" Struct used to define a 2D vector for Paraview.

The y-component of velocity is flipped for compatibility with Paraview.
Also, a zero magnitude z-component of velocity is added so that Paraview
will recognize the array as a vector.
"""
struct Vector2D
    name::String
    jld_dataname_vxyz::String
    jld_dataname_vx::String
    jld_dataname_vy::String
    jld_dataname_vmag::String
    number_type::String
    units::String
    header::String
    grid_type::String
    vectorx::Matrix{Float64}
    vectory::Matrix{Float64}
    vectorxyz::Array{Float64,3}
    vmag::Matrix{Float64}
end

function Vector2D(
    name::String, 
    jld_dataname_vxyz::String, 
    jld_dataname_vx::String,
    jld_dataname_vy::String, 
    jld_dataname_vmag::String, 
    number_type::String,
    units::String, 
    header::String, 
    grid_type::String,
    vectorx::Matrix{Float64}, 
    vectory::Matrix{Float64}
)
    vectorxyz = make_vector(vectorx, vectory)
    vmag = calculate_velocity_magnitude(
        size(vectorx, 1), size(vectorx, 2), vectorx, vectory)
    return Vector2D(
        name, jld_dataname_vxyz, jld_dataname_vx, jld_dataname_vy,
        jld_dataname_vmag, number_type, units, header, grid_type,
        vectorx, vectory, vectorxyz, vmag
    )
end

""" Make a 3-component array using array stacking.
"""
function make_vector(
    vectorx::Matrix{Float64}, 
    vectory::Matrix{Float64}
)::Array{Float64,3}
    ny, nx = size(vectorx)
    vectorz = zeros(Float64, ny, nx)
    vectorxyz = cat(vectorx, -vectory, vectorz, dims=3)
    return vectorxyz
end

""" Calculate the magnitude of velocity vectors.
"""
function calculate_velocity_magnitude(
    ynum::Int, 
    xnum::Int,
    vectorx::Matrix{Float64},
    vectory::Matrix{Float64}
)::Matrix{Float64}
    vmag = zeros(Float64, ynum, xnum)
    for i in 1:ynum
        for j in 1:xnum
            vyy = vectory[i, j]
            vxx = vectorx[i, j]
            vmag[i, j] = sqrt(vxx*vxx + vyy*vyy)
        end
    end
    return vmag
end


"""
Struct used to store jld2 information associated with markers.
"""
struct Markers2djld
    jld_markerfile::String
    jld_dataname_xy::String
    jld_dataname_ids::String
    nmarkers::Int
    marker_xy_km_array::Matrix{Float64}
    marker_id_array::Vector{Float64}
end

function Markers2djld(marker_data::Dict{String, Any})
    return Markers2djld(
        marker_data["jld_markerfile"],
        marker_data["jld_dataname_xy"],
        marker_data["jld_dataname_ids"],
        marker_data["nmarkers"],
        marker_data["marker_xy_km_array"],
        marker_data["marker_id_array"]
    )
end

""" Struct used to store jld2 information associated with 2D topography.
"""
struct Topo2djld
    jld_markerfile::String
    jld_dataname_xy::String
    jld_dataname_ids::String
    ntopo::Int
    topo_xy_km_array::Matrix{Float64}
    topo_id_array::Vector{Float64}
end

function Topo2djld(topo_data::Dict{String, Any})
    return Topo2djld(
        topo_data["jld_markerfile"],
        topo_data["jld_dataname_xy"],
        topo_data["jld_dataname_ids"],
        topo_data["ntopo"],
        topo_data["topo_xy_km_array"],
        topo_data["topo_id_array"]
    )
end

end # module 