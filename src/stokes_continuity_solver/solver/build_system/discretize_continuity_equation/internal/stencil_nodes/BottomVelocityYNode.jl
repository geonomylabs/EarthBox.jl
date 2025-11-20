module BottomVelocityYNode

import EarthBox.Domain: basic_cell_along_lower_boundary
import EarthBox.Domain: basic_cell_along_upper_boundary
import ....BuildStructs: StokesBuildData
import ....BuildStructs: CellIndices

function coefficients_and_rhs_term(
    inz::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Int
    i = cell_indices.i
    j = cell_indices.j
    ipr = cell_indices.ipr
    ynum = build_data.grid.ynum
    ystp = build_data.grid.ystp
    btopy = build_data.bc.btopy
    pscale = build_data.pscale
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    if !basic_cell_along_lower_boundary(i, ynum)
        ivyD = ipr - 1
        Lv[inz] = pscale / ystp[i]
        Li[inz] = ipr
        Lj[inz] = ivyD
        Lii[inz] = i
        Ljj[inz] = j
        inz_yD = inz
        inz = inz + 1
        if basic_cell_along_upper_boundary(i)
            Lv[inz_yD] = Lv[inz_yD] - btopy[j+1, 2] * pscale / ystp[i]
            R[ipr] = R[ipr] + btopy[j+1, 1] * pscale / ystp[i]
        end
    end
    return inz
end

end # module 