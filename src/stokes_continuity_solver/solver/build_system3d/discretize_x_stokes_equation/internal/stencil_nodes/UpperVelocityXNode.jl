module UpperVelocityXNode

import EarthBox.Domain: basic_cell_along_upper_boundary
import ....BuildStructs: StokesBuildData
import ....BuildStructs: CellIndices

function coefficients_and_rhs_term(
    inz::Int,
    inz_c::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Int
    ivx = cell_indices.ivx
    i = cell_indices.i
    j = cell_indices.j
    k = cell_indices.k
    ystp = build_data.grid.ystp
    ystpc = build_data.grid.ystpc
    etas_xy = build_data.etas_xy
    btopx = build_data.bc.btopx
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    if !basic_cell_along_upper_boundary(i)
        # The within-cell stride between vx of cell (i,j,k) and vx of cell
        # (i-1,j,k) is 4 in 3D.
        ivxU = ivx - 4
        Lv[inz] = etas_xy[i, j+1, k] / ystpc[i] / ystp[i]
        Li[inz] = ivx
        Lj[inz] = ivxU
        Lii[inz] = i
        Ljj[inz] = j
        inz = inz + 1
    else
        # vxU is an upper ghost boundary node.
        Lv[inz_c] = (
            Lv[inz_c] +
            btopx[j+1, k+1, 2] * etas_xy[i, j+1, k] / ystpc[i] / ystp[i]
        )
        R[ivx] = (
            R[ivx] -
            btopx[j+1, k+1, 1] * etas_xy[i, j+1, k] / ystpc[i] / ystp[i]
        )
    end
    return inz
end

end # module
