""" Functions used to interpolate topography array info at marker locations.

Interpolated information:
- Topography elevation: gridt[2, :] (meters)
- Extrusion thickness: gridt[7, :] (meters)
"""
module TopoInterpolation

import EarthBox.ModelDataContainer: ModelData

"""
    find_closest_topography_node_to_marker(x_marker, gridt, dx_topo, toponum)

Find index of closest topography node to marker.

# Arguments
- `x_marker::Float64`: X-coordinate of marker location in meters.
- `gridt::Array{Float64,2}`: Multi-dimensional topography grid array with elevation, 
  velocity, antidiffusion information and extrusion thickness associated with volcanic 
  processes. In this function the x-coordinates stored in gridt[1, :] are used to find 
  the closest node to the marker.
- `dx_topo::Float64`: Grid spacing of topography nodes.
- `toponum::Int`: Number of topography nodes.

# Returns
- `ixn::Int`: Index of closest topography node
- `dx_dist::Float64`: Distance to closest topography node.
"""
@inline function find_closest_topography_node_to_marker(
    x_marker::Float64,
    gridt::Array{Float64,2},
    dx_topo::Float64,
    toponum::Int
)::Tuple{Int,Float64}
    ixn = floor(Int, (x_marker - gridt[1, 1])/dx_topo - 0.5)
    ixn = max(ixn, 1)
    if ixn > toponum - 1
        ixn = toponum - 1
    end
    dx_dist = (x_marker - gridt[1, ixn])/dx_topo
    return ixn, dx_dist
end

"""
    calculate_topography(dx_dist, gridt, ixn)

Compute topography elevation at marker location.

# Arguments
- `dx_dist::Float64`: Distance to closest topography node.
- `gridt::Array{Float64,2}`: Multi-dimensional topography grid array with elevation, 
  velocity, antidiffusion information and extrusion thickness associated with volcanic 
  processes. In this function the elevation stored in gridt[2, :] is used to interpolate 
  at off-node locations.
- `ixn::Int`: Index of closest topography node.

# Returns
- `y_topo::Float64`: Y-coordinate of topography elevation in meters.
"""
@inline function calculate_topography(
    dx_dist::Float64,
    gridt::Array{Float64,2},
    ixn::Int
)::Float64
    y_topo = gridt[2, ixn]*(1.0 - dx_dist) + gridt[2, ixn+1]*dx_dist
    return y_topo
end

"""
    calculate_extrusion_thickness(dx_dist, gridt, ixn)

Compute extrusion thickness at marker location.

# Arguments
- `dx_dist::Float64`: Distance to closest topography node.
- `gridt::Array{Float64,2}`: Multi-dimensional topography grid array with elevation, 
  velocity, antidiffusion information and extrusion thickness associated with volcanic 
  processes. In this function the extrusion thickness stored in gridt[7, :] is used to 
  interpolate at off-node locations.
- `ixn::Int`: Index of closest topography node.

# Returns
- `ext_thick::Float64`: Extrusion thickness at marker location in meters.
"""
@inline function calculate_extrusion_thickness(
    dx_dist::Float64,
    gridt::Array{Float64,2},
    ixn::Int
)::Float64
    ext_thick = gridt[7, ixn]*(1.0 - dx_dist) + gridt[7, ixn+1]*dx_dist
    return ext_thick
end 

end