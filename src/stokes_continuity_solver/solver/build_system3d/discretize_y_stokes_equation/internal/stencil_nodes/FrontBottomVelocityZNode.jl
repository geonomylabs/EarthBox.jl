module FrontBottomVelocityZNode

import ...Predicates3D: basic_cell_along_front_boundary
import ....BuildStructs: StokesBuildData
import ....BuildStructs: CellIndices

function coefficients_and_rhs_term(
    inz::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData,
    drhodz_gy_dt::Float64
)::Tuple{Int,Int,Int}
    i = cell_indices.i
    j = cell_indices.j
    k = cell_indices.k
    ivy = cell_indices.ivy
    dshift_to_vxF = build_data.grid.dshift_to_vxF
    ystpc = build_data.grid.ystpc
    zstp = build_data.grid.zstp
    etas_yz = build_data.etas_yz
    bfrontz = build_data.bc.bfrontz
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    inz_FD = -1
    inz_BD = -1
    if !basic_cell_along_front_boundary(k)
        # vz of cell (i+1, j, k-1). Cross-i stride is 4; cross-k is dshift.
        ivzFD = ivy + 5 - dshift_to_vxF
        Lv[inz] = -etas_yz[i+1, j, k] / ystpc[i+1] / zstp[k] - drhodz_gy_dt / 4.0
        Li[inz] = ivy
        Lj[inz] = ivzFD
        Lii[inz] = i
        Ljj[inz] = j
        inz_FD = inz
        inz = inz + 1
    else
        # vzFD is at the front boundary. Substitute via bfrontz; vzBD lives
        # in cell (i+1, j, k).
        ivzBD = ivy + 5
        Lv[inz] = bfrontz[i+2, j+1, 2] * (
            -etas_yz[i+1, j, k] / ystpc[i+1] / zstp[k] - drhodz_gy_dt / 4.0
        )
        Li[inz] = ivy
        Lj[inz] = ivzBD
        Lii[inz] = i
        Ljj[inz] = j
        inz_BD = inz
        inz = inz + 1
        inz_FD = 0
        R[ivy] = R[ivy] - bfrontz[i+2, j+1, 1] * (
            -etas_yz[i+1, j, k] / ystpc[i+1] / zstp[k] - drhodz_gy_dt / 4.0
        )
    end
    return inz, inz_FD, inz_BD
end

end # module
