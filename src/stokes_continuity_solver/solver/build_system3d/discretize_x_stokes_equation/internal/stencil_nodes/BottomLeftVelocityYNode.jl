module BottomLeftVelocityYNode

import EarthBox.Domain: basic_cell_along_upper_boundary
import EarthBox.Domain: basic_cell_along_lower_boundary
import ....BuildStructs: StokesBuildData
import ....BuildStructs: CellIndices

function coefficients_and_rhs_term(
    inz::Int,
    inz_UL::Int,
    inz_LL::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Int
    ivx = cell_indices.ivx
    i = cell_indices.i
    j = cell_indices.j
    k = cell_indices.k
    ynum = build_data.grid.ynum
    xstpc = build_data.grid.xstpc
    ystp = build_data.grid.ystp
    etas_xy = build_data.etas_xy
    bbottomy = build_data.bc.bbottomy
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    if !basic_cell_along_lower_boundary(i, ynum)
        ivyLL = ivx + 1
        if !basic_cell_along_upper_boundary(i)
            Lv[inz] = -etas_xy[i+1, j+1, k] / xstpc[j+1] / ystp[i]
            Li[inz] = ivx
            Lj[inz] = ivyLL
            Lii[inz] = i
            Ljj[inz] = j
            inz = inz + 1
        else
            # Basic node is at top boundary and a coefficient was already
            # defined for vyLL by UpperLeftVelocityYNode. Update it.
            Lv[inz_LL] = (
                Lv[inz_LL] -
                etas_xy[i+1, j+1, k] / xstpc[j+1] / ystp[i]
            )
        end
    else
        # vyLL is a ghost boundary node. Substitute the boundary equation
        # that relates vyLL to vyUL, then update the coefficient on vyUL
        # and the rhs term of the central node.
        Lv[inz_UL] = (
            Lv[inz_UL] -
            bbottomy[j+1, k+1, 2] * etas_xy[i+1, j+1, k] / xstpc[j+1] / ystp[i]
        )
        R[ivx] = (
            R[ivx] +
            bbottomy[j+1, k+1, 1] * etas_xy[i+1, j+1, k] / xstpc[j+1] / ystp[i]
        )
    end
    return inz
end

end # module
