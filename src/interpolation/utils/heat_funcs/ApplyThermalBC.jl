module ApplyThermalBC

"""
    apply_thermal_bc_along_boundaries!(model)

Apply thermal boundary conditions along boundaries of the basic grid.

This function is used after interpolating marker information to the basic
grid to ensure that boundary condition values are correct along boundary
nodes.

# Arguments
- `model`: ModelData struct containing grid and boundary condition information

# Updated Arrays
## Updated arrays from group `model.heat_equation.arrays.temperature`
- `tk1.array::Matrix{Float64}`:
    - Marker temperature (K) interpolated to basic grid
"""
function apply_thermal_bc_along_boundaries!(model)
    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value
    btopt = model.bcs.arrays.temperature.btopt.array
    bbottomt = model.bcs.arrays.temperature.bbottomt.array
    bleftt = model.bcs.arrays.temperature.bleftt.array
    brightt = model.bcs.arrays.temperature.brightt.array
    tk1 = model.heat_equation.arrays.temperature.tk1.array

    # Upper and lower boundaries
    for j in 2:xnum-1
        # Upper boundary
        tk1[1, j] = btopt[j, 1] + btopt[j, 2] * tk1[2, j]
        # Lower boundary
        tk1[ynum, j] = bbottomt[j, 1] + bbottomt[j, 2] * tk1[ynum-1, j]
    end

    # Left and right boundaries
    for i in 1:ynum
        # Left boundary
        tk1[i, 1] = bleftt[i, 1] + bleftt[i, 2] * tk1[i, 2]
        # Right boundary
        tk1[i, xnum] = brightt[i, 1] + brightt[i, 2] * tk1[i, xnum-1]
    end
end

end # module 