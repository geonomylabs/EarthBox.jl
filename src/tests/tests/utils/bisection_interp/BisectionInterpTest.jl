module BisectionInterpTest

import EarthBox.MathTools: linear_interp_bisection

function run_test()::Nothing
    gridx = collect(range(0.0, 500_000.0, length=100))
    nx = length(gridx)
    gridy = collect(range(3.0, 10.0, length=nx))

    # print min and max y
    println("min y = ", minimum(gridy))
    println("max y = ", maximum(gridy))

    x = 800_000.0
    y = linear_interp_bisection(gridx, gridy, x)
    println("x ", x, " -> y = ", y)

    return nothing
end

end # module 