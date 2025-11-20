module ConstantGradientBC

import EarthBox.ModelDataContainer: ModelData

"""
    set_constant_vertical_gradient_at_top_boundary!(
        ystp::Float64,
        dtdy::Float64,
        model::ModelData
    )::Nothing

Define coefficients for constant gradient along the top boundary.

# Arguments
- `ystp::Float64`: Vertical step size
- `dtdy::Float64`: Temperature gradient in y-direction
- `model::ModelData`: Model data structure

# Updated array from group bcs.arrays.temperature
- btopt.array::Matrix{Float64}

# Background
The following equation is used to define temperature at basic nodes along the 
top boundary:

tk[1][j] = btopt[j][0] + tk[2][j]*btop[j][1]

where j is the horizontal index of the basic grid.

Temperature at the top is set equal to the temperature of the next deepest 
basic grid node less the temperature change assuming a constant gradient (dtdy) 
across the basic grid cell with thickness ystp:

tk[1][j] = -ystp*dtdy + tk[2][j]

where btopt[j][1] = -ystp*dtdy and btopt[j][2] = 1.
"""
function set_constant_vertical_gradient_at_top_boundary!(
    ystp::Float64,
    dtdy::Float64,
    model::ModelData
)::Nothing
    xnum = model.grids.parameters.geometry.xnum.value
    btopt = model.bcs.arrays.temperature.btopt.array
    for j in 1:xnum
        btopt[j, 1] = -ystp * dtdy
        btopt[j, 2] = 1.0
    end
    return nothing
end

"""
    set_constant_vertical_gradient_at_bottom_boundary!(
        ystp::Float64,
        dtdy::Float64,
        model::ModelData
    )::Nothing

Define coefficients for constant gradient along the bottom boundary.

# Arguments
- `ystp::Float64`: Vertical step size
- `dtdy::Float64`: Temperature gradient in y-direction
- `model::ModelData`: Model data structure

# Updated array from group bcs.arrays.temperature
- bbottomt.array::Matrix{Float64}

# Background
The following equation is used to define temperature at basic nodes along the 
bottom boundary:

tk[ynum][j] = bbottomt[j][0] + tk[ynum-1][j]*bbottomt[j][1]

where j is the horizontal index of the basic grid.

Temperature at the bottom is set equal to the temperature of the next shallowest 
basic grid node plus the temperature change assuming a constant gradient (dtdy) 
across the basic grid cell with thickness ystp:

tk[ynum][j] = ystp*dtdy + tk[ynum-1][j]

where bbottomt[j][1] = ystp*dtdy and bbottomt[j][2] = 1.
"""
function set_constant_vertical_gradient_at_bottom_boundary!(
    ystp::Float64,
    dtdy::Float64,
    model::ModelData
)::Nothing
    xnum = model.grids.parameters.geometry.xnum.value
    bbottomt = model.bcs.arrays.temperature.bbottomt.array
    for j in 1:xnum
        bbottomt[j, 1] = ystp * dtdy
        bbottomt[j, 2] = 1.0
    end
    return nothing
end

end # module 