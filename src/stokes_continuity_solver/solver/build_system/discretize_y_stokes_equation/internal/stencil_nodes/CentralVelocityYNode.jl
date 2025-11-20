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
    ivy = cell_indices.ivy
    xstp = build_data.grid.xstp
    ystp = build_data.grid.ystp
    xstpc = build_data.grid.xstpc
    ystpc = build_data.grid.ystpc
    etan = build_data.etan
    etas = build_data.etas
    RY = build_data.rhs.RY
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    Lv[inz] = (
        -2*(etan[i+1,j]/ystp[i+1] + etan[i,j]/ystp[i])/ystpc[i+1] 
        - (etas[i+1,j+1]/xstpc[j+1] + etas[i+1,j]/xstpc[j])/xstp[j] 
        - drhody_gy_dt
    )
    Li[inz] = ivy
    Lj[inz] = ivy
    Lii[inz] = i
    Ljj[inz] = j
    inz_c = inz
    inz = inz + 1
    R[ivy] = RY[i+1,j+1]
    return inz, inz_c
end

end # module 