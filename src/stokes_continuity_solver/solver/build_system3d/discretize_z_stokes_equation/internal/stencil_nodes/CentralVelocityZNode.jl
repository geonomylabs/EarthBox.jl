module CentralVelocityZNode

import ....BuildStructs: StokesBuildData
import ....BuildStructs: CellIndices

function coefficients_and_rhs_term(
    inz::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Tuple{Int,Int}
    i = cell_indices.i
    j = cell_indices.j
    k = cell_indices.k
    ivz = cell_indices.ivz
    xstp = build_data.grid.xstp
    ystp = build_data.grid.ystp
    zstp = build_data.grid.zstp
    xstpc = build_data.grid.xstpc
    ystpc = build_data.grid.ystpc
    zstpc = build_data.grid.zstpc
    etan = build_data.etan
    etas_xz = build_data.etas_xz
    etas_yz = build_data.etas_yz
    RZ = build_data.rhs.RZ
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    Lv[inz] = (
        -2.0 * (etan[i, j, k+1] / zstp[k+1] + etan[i, j, k] / zstp[k]) / zstpc[k+1]
        - (etas_xz[i, j+1, k+1] / xstpc[j+1] + etas_xz[i, j, k+1] / xstpc[j]) / xstp[j]
        - (etas_yz[i+1, j, k+1] / ystpc[i+1] + etas_yz[i, j, k+1] / ystpc[i]) / ystp[i]
    )
    Li[inz] = ivz
    Lj[inz] = ivz
    Lii[inz] = i
    Ljj[inz] = j
    inz_c = inz
    inz = inz + 1
    R[ivz] = RZ[i+1, j+1, k+1]
    return inz, inz_c
end

end # module
