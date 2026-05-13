module BackUpperVelocityZNode

import ...Predicates3D: basic_cell_along_front_boundary
import ...Predicates3D: basic_cell_along_back_boundary
import ....BuildStructs: StokesBuildData
import ....BuildStructs: CellIndices

function coefficients_and_rhs_term(
    inz::Int,
    inz_FU::Int,
    inz_BU::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData,
    drhodz_gy_dt::Float64
)::Int
    i = cell_indices.i
    j = cell_indices.j
    k = cell_indices.k
    ivy = cell_indices.ivy
    znum = build_data.grid.znum
    dshift_to_vxF = build_data.grid.dshift_to_vxF
    ystpc = build_data.grid.ystpc
    zstp = build_data.grid.zstp
    etas_yz = build_data.etas_yz
    bbackz = build_data.bc.bbackz
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    if !basic_cell_along_back_boundary(k, znum)
        # vz of cell (i, j, k) — same cell as vyC.
        ivzBU = ivy + 1
        if !basic_cell_along_front_boundary(k)
            Lv[inz] = -etas_yz[i+1, j, k+1] / ystpc[i+1] / zstp[k] - drhodz_gy_dt / 4.0
            Li[inz] = ivy
            Lj[inz] = ivzBU
            Lii[inz] = i
            Ljj[inz] = j
            inz = inz + 1
        else
            # Basic node is at front boundary; vzBU coefficient was already
            # defined by FrontUpperVelocityZNode at inz_BU. Update it.
            Lv[inz_BU] = (
                Lv[inz_BU] -
                etas_yz[i+1, j, k+1] / ystpc[i+1] / zstp[k] - drhodz_gy_dt / 4.0
            )
        end
    else
        # vzBU is a back ghost boundary node. Substitute via bbackz.
        Lv[inz_FU] = Lv[inz_FU] + bbackz[i+1, j+1, 2] * (
            -etas_yz[i+1, j, k+1] / ystpc[i+1] / zstp[k] - drhodz_gy_dt / 4.0
        )
        R[ivy] = R[ivy] - bbackz[i+1, j+1, 1] * (
            -etas_yz[i+1, j, k+1] / ystpc[i+1] / zstp[k] - drhodz_gy_dt / 4.0
        )
    end
    return inz
end

end # module
