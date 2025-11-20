module MarkerDataContainer

mutable struct MarkerData
    # Number of markers
    marknum::Union{Int, Nothing}
    # Marker material id's
    marker_matid::Union{Vector{Int16}, Nothing}
    # X-coordinates of markers in meters
    marker_x::Union{Vector{Float64}, Nothing}
    # Y-coordinates of markers in meters
    marker_y::Union{Vector{Float64}, Nothing}
end

# Constructor with default values
MarkerData(; marknum=nothing, marker_matid=nothing, marker_x=nothing,
    marker_y=nothing) = MarkerData(marknum, marker_matid, marker_x, marker_y)

end # module 