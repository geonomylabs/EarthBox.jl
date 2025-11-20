module ZeroHeatFluxBC

using EarthBox.ModelDataContainer: ModelData

"""
    set_zero_heat_flux_at_left_boundary!(model::ModelData)::Nothing

Define coefficients for zero flux along left boundary.

# Arguments
- `model::ModelData`: The model data structure

# Updated array from group bcs.arrays.temperature
- bleftt.array::Matrix{Float64}

# Background
The following equation is used to define temperature at basic nodes along the 
left boundary:

tk[i][1] = bleftt[i][1] + bleftt[i][2]*tk[i][2]

where i is the vertical index of the basic grid.

Zero flux is implemented by setting temperature along the boundary equal to
the next internal grid node toward the right.
"""
function set_zero_heat_flux_at_left_boundary!(model::ModelData)::Nothing
    ynum = model.grids.parameters.geometry.ynum.value
    bleftt = model.bcs.arrays.temperature.bleftt.array
    for i in 1:ynum
        bleftt[i, 1] = 0.0
        bleftt[i, 2] = 1.0
    end
    return nothing
end

"""
    set_zero_heat_flux_at_right_boundary!(model::ModelData)::Nothing

Define coefficients for zero flux along right boundary.

# Arguments
- `model::ModelData`: The model data structure

# Updated array from group bcs.arrays.temperature
- brightt.array::Matrix{Float64}

# Background
The following equation is used to define temperature:

tk[i][xnum] = brightt[i][1] + brightt[i][2]*tk[i][xnum-1]

where i is the vertical index of the basic grid.

Zero flux is implemented by setting temperature along the boundary equal to
the next internal grid node toward the left.
"""
function set_zero_heat_flux_at_right_boundary!(model::ModelData)::Nothing
    ynum = model.grids.parameters.geometry.ynum.value
    brightt = model.bcs.arrays.temperature.brightt.array
    for i in 1:ynum
        brightt[i, 1] = 0.0
        brightt[i, 2] = 1.0
    end
    return nothing
end

end # module