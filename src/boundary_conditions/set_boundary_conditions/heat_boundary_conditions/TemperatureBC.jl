module TemperatureBC

using EarthBox.ModelDataContainer: ModelData

"""
    set_temperature_at_top_boundary!(
        temperature_kelvins::Float64, model::ModelData)::Nothing

Define coefficients for prescribed temperature along the top boundary.

# Arguments
- `temperature_kelvins::Float64`: The temperature in Kelvin to set at the top boundary
- `model::ModelData`: The model data structure

# Background
The following equation is used to define temperature at basic nodes along the 
top boundary:

tk[1][j] = btopt[j][1] + tk[2][j]*btop[j][2]

where j is the horizontal index of the basic grid.

Temperature is set to a constant value at the top boundary by setting
btopt[j][1] to the prescribed temperature and btopt[j][2] to zero:

tk[1][j] = btopt[j][1]

"""
function set_temperature_at_top_boundary!(
    temperature_kelvins::Float64, 
    model::ModelData
)::Nothing
    xnum = model.grids.parameters.geometry.xnum.value
    btopt = model.bcs.arrays.temperature.btopt.array
    for j in 1:xnum
        btopt[j, 1] = temperature_kelvins
        btopt[j, 2] = 0.0
    end
    return nothing
end

"""
    set_temperature_at_bottom_boundary!(
        temperature_kelvins::Float64, model::ModelData)::Nothing

Define coefficients for prescribed temperature along the bottom boundary.

# Arguments
- `temperature_kelvins`: The temperature in Kelvin to set at the bottom boundary
- `model`: The model data structure

# Updated array from group bcs.arrays.temperature
- bbottomt.array::Matrix{Float64}

# Background
The following equation is used to define temperature at basic nodes along the 
bottom boundary:

tk[ynum][j] = bbottomt[j][1] + tk[ynum-1][j]*bbottomt[j][2]

where j is the horizontal index of the basic grid.

Temperature is set to a constant value at the bottom boundary by setting
bbottomt[j][1] to the prescribed temperature and bbottomt[j][2] to zero:

tk[ynum][j] = bbottomt[j][1]

"""
function set_temperature_at_bottom_boundary!(
    temperature_kelvins::Float64, 
    model::ModelData
)::Nothing
    xnum = model.grids.parameters.geometry.xnum.value
    bbottomt = model.bcs.arrays.temperature.bbottomt.array
    for j in 1:xnum
        bbottomt[j, 1] = temperature_kelvins
        bbottomt[j, 2] = 0.0
    end
    return nothing
end

"""
    set_temperature_at_left_boundary(temperature_kelvins::Float64, model::ModelData)

Define coefficients for prescribed temperature along the left boundary.

# Arguments
- `temperature_kelvins`: The temperature in Kelvin to set at the left boundary
- `model`: The model data structure

# Updated array from group bcs.arrays.temperature
- bleftt.array::Matrix{Float64}

# Background
The following equation is used to define temperature at basic nodes along the 
left boundary:

tk[i][1] = bleftt[i][1] + bleftt[i][2]*tk[i][2]

where i is the vertical index of the basic grid.

Temperature is set to a constant value at the left boundary by setting
bleftt[i][0] to the prescribed temperature and bleftt[i][1] to zero:

tk[i][1] = bleftt[i][1]

"""
function set_temperature_at_left_boundary!(temperature_kelvins::Float64, model::ModelData)::Nothing
    ynum = model.grids.parameters.geometry.ynum.value
    bleftt = model.bcs.arrays.temperature.bleftt.array
    for i in 1:ynum
        bleftt[i, 1] = temperature_kelvins
        bleftt[i, 2] = 0.0
    end
    return nothing
end

"""
    set_variable_temperature_at_left_boundary!(
        temperature_top_kelvins::Float64, 
        ystp::Float64, 
        dtdy::Float64, 
        model::ModelData
    )::Nothing

Define coefficients for variable vertical temperature at left boundary.

# Arguments
- `temperature_top_kelvins`: The temperature in Kelvin at the top of the left boundary
- `ystp`: The average vertical grid spacing
- `dtdy`: The temperature gradient along the y-axis
- `model`: The model data structure

# Updated array from group bcs.arrays.temperature
- bleftt.array::Matrix{Float64}

# Background
The following equation is used to define temperature at basic nodes along the 
left boundary:

tk[i][1] = bleftt[i][1] + bleftt[i][2]*tk[i][2]

where i is the vertical index of the basic grid.

Temperature is set to a variable value at the left boundary by setting
bleftt[i][1] to the variable temperature and bleftt[i][2] to zero:

tk[i][1] = bleftt[i][1] + float(i)*ystp*dtdy
"""
function set_variable_temperature_at_left_boundary!(
    temperature_top_kelvins::Float64, 
    ystp::Float64, 
    dtdy::Float64, 
    model::ModelData
)::Nothing
    ynum = model.grids.parameters.geometry.ynum.value
    bleftt = model.bcs.arrays.temperature.bleftt.array
    for i in 1:ynum
        bleftt[i, 1] = temperature_top_kelvins + (i-1)*ystp*dtdy
        bleftt[i, 2] = 0.0
    end
    return nothing
end

"""
    set_temperature_at_right_boundary!(
        temperature_kelvins::Float64, model::ModelData)::Nothing

Define coefficients for prescribed temperature along the right boundary.

# Arguments
- `temperature_kelvins`: The temperature in Kelvin to set at the right boundary
- `model`: The model data structure

# Updated array from group bcs.arrays.temperature
- brightt.array::Matrix{Float64}

# Background
The following equation is used to define temperature at basic nodes along the 
right boundary:

tk[i][xnum] = brightt[i][1] + brightt[i][2]*tk[i][xnum-1]

where i is the vertical index of the basic grid.

Temperature is set to a constant value at the right boundary by setting
brightt[i][1] to the prescribed temperature and brightt[i][2] to zero:

tk[i][xnum] = brightt[i][1]
"""
function set_temperature_at_right_boundary!(
    temperature_kelvins::Float64,
    model::ModelData
)::Nothing
    ynum = model.grids.parameters.geometry.ynum.value
    brightt = model.bcs.arrays.temperature.brightt.array
    for i in 1:ynum
        brightt[i, 1] = temperature_kelvins
        brightt[i, 2] = 0.0
    end
    return nothing
end

"""
    set_variable_temperature_at_right_boundary!(
        temperature_top::Float64, 
        ystp::Float64, 
        dtdy::Float64, 
        model::ModelData
    )::Nothing

Define coefficients for variable temperature along the right boundary.

# Arguments
- `temperature_top`: The temperature at the top of the right boundary
- `ystp`: The average vertical grid spacing
- `dtdy`: The temperature gradient along the y-axis
- `model`: The model data structure

# Updated array from group bcs.arrays.temperature
- brightt.array::Matrix{Float64}

# Background
The following equation is used to define temperature at basic nodes along the 
right boundary:

tk[i][xnum] = brightt[i][1] + brightt[i][2]*tk[i][xnum-1]

where i is the vertical index of the basic grid.

Temperature is set to a variable value at the right boundary by setting
brightt[i][1] to the variable temperature and brightt[i][2] to zero:

tk[i][xnum] = brightt[i][1] + float(i)*ystp*dtdy
"""
function set_variable_temperature_at_right_boundary!(
    temperature_top::Float64, 
    ystp::Float64, 
    dtdy::Float64, 
    model::ModelData
)::Nothing
    ynum = model.grids.parameters.geometry.ynum.value
    brightt = model.bcs.arrays.temperature.brightt.array
    for i in 1:ynum
        brightt[i, 1] = temperature_top + (i-1)*ystp*dtdy
        brightt[i, 2] = 0.0
    end
    return nothing
end

end # module