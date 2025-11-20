module PoissonCoefficientsManager

import ...GridDataManager: GridData

struct PoissonCoefficients
    xkf::Float64
    ykf::Float64
    zkf::Float64
    xkf2::Float64
    ykf2::Float64
    zkf2::Float64
    xykf::Float64
    xzkf::Float64
    yzkf::Float64
end

function PoissonCoefficients(
    n::Int64,
    grid_data::GridData,
)::PoissonCoefficients
    xstp = grid_data.xstp[n]
    ystp = grid_data.ystp[n]
    zstp = grid_data.zstp[n]

    xkf = 1.0 / xstp^2
    ykf = 1.0 / ystp^2
    zkf = 1.0 / zstp^2
    xkf2 = 2.0 / xstp^2
    ykf2 = 2.0 / ystp^2
    zkf2 = 2.0 / zstp^2
    xykf = 1.0 / xstp / ystp
    xzkf = 1.0 / xstp / zstp
    yzkf = 1.0 / ystp / zstp

    return PoissonCoefficients(xkf, ykf, zkf, xkf2, ykf2, zkf2, xykf, xzkf, yzkf)
end

end # module