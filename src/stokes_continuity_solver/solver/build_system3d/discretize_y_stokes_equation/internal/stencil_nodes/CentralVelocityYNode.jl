module CentralVelocityYNode

import ....BuildStructs: StokesBuildData
import ....BuildStructs: CellIndices

function coefficients_and_rhs_term(
    inz::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData,
    drhody_gy_dt::Float64
)::Tuple{Int,Int}
    i = cell_indices.i
    j = cell_indices.j
    k = cell_indices.k
    ivy = cell_indices.ivy
    xstp = build_data.grid.xstp
    ystp = build_data.grid.ystp
    zstp = build_data.grid.zstp
    xstpc = build_data.grid.xstpc
    ystpc = build_data.grid.ystpc
    zstpc = build_data.grid.zstpc
    etan = build_data.etan
    etas_xy = build_data.etas_xy
    etas_yz = build_data.etas_yz
    RY = build_data.rhs.RY
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    Lv[inz] = (
        -2.0 * (etan[i+1, j, k] / ystp[i+1] + etan[i, j, k] / ystp[i]) / ystpc[i+1]
        - (etas_xy[i+1, j+1, k] / xstpc[j+1] + etas_xy[i+1, j, k] / xstpc[j]) / xstp[j]
        - (etas_yz[i+1, j, k+1] / zstpc[k+1] + etas_yz[i+1, j, k] / zstpc[k]) / zstp[k]
        - drhody_gy_dt
    )
    Li[inz] = ivy
    Lj[inz] = ivy
    Lii[inz] = i
    Ljj[inz] = j
    inz_c = inz
    inz = inz + 1
    R[ivy] = RY[i+1, j+1, k+1]
    return inz, inz_c
end

end # module
