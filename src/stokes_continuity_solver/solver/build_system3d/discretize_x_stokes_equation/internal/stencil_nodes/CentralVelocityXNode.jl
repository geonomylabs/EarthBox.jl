module CentralVelocityXNode

import ....BuildStructs: StokesBuildData
import ....BuildStructs: CellIndices

function coefficients_and_rhs_term(
    inz::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Tuple{Int,Int}
    ivx = cell_indices.ivx
    i = cell_indices.i
    j = cell_indices.j
    k = cell_indices.k
    xstp = build_data.grid.xstp
    ystp = build_data.grid.ystp
    zstp = build_data.grid.zstp
    xstpc = build_data.grid.xstpc
    ystpc = build_data.grid.ystpc
    zstpc = build_data.grid.zstpc
    etan = build_data.etan
    etas_xy = build_data.etas_xy
    etas_xz = build_data.etas_xz
    RX = build_data.rhs.RX
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    Lv[inz] = (
        -2.0 * (etan[i, j+1, k] / xstp[j+1] + etan[i, j, k] / xstp[j]) / xstpc[j+1]
        - (etas_xy[i+1, j+1, k] / ystpc[i+1] + etas_xy[i, j+1, k] / ystpc[i]) / ystp[i]
        - (etas_xz[i, j+1, k+1] / zstpc[k+1] + etas_xz[i, j+1, k] / zstpc[k]) / zstp[k]
    )
    Li[inz] = ivx
    Lj[inz] = ivx
    Lii[inz] = i
    Ljj[inz] = j
    # Store the non-zero index of the central node so that vxC coefficients
    # and rhs term can be updated when stencil includes ghost boundary nodes.
    inz_c = inz
    inz = inz + 1
    R[ivx] = RX[i+1, j+1, k+1]
    return inz, inz_c
end

end # module
