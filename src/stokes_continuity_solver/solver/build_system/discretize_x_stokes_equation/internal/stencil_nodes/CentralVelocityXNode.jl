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
    xstp = build_data.grid.xstp
    ystp = build_data.grid.ystp
    xstpc = build_data.grid.xstpc
    ystpc = build_data.grid.ystpc
    etan = build_data.etan
    etas = build_data.etas
    RX = build_data.rhs.RX
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    Lv[inz] = (
        -2.0 * (etan[i, j+1] / xstp[j+1] + etan[i, j] / xstp[j]) / xstpc[j+1] -
        (etas[i+1, j+1] / ystpc[i+1] + etas[i, j+1] / ystpc[i]) / ystp[i]
    )
    Li[inz] = ivx
    Lj[inz] = ivx
    Lii[inz] = i
    Ljj[inz] = j
    # Store the non-zero index of the central node so that vxC coefficients and
    # rhs term can be updated when stencil includes ghost boundary nodes.
    inz_c = inz
    inz = inz + 1
    # Define the right-hand side term.
    R[ivx] = RX[i+1, j+1]
    return inz, inz_c
end

end # module 