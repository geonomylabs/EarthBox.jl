module BottomVelocityXNode

import EarthBox.Domain: basic_cell_along_lower_boundary
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
    ynum = build_data.grid.ynum
    ystp = build_data.grid.ystp
    ystpc = build_data.grid.ystpc
    etas = build_data.etas
    bbottomx = build_data.bc.bbottomx
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    if !basic_cell_along_lower_boundary(i, ynum)
        ivxD = ivx + 3
        Lv[inz] = etas[i+1, j+1] / ystpc[i+1] / ystp[i]
        Li[inz] = ivx
        Lj[inz] = ivxD
        Lii[inz] = i
        Ljj[inz] = j
        inz = inz + 1
    else
        # Apply boundary conditions since vxD is a ghost boundary node
        # Adding new term to vxC coefficient
        Lv[inz_c] = (
            Lv[inz_c] +
            bbottomx[j+1, 2] * etas[i+1, j+1] / ystpc[i+1] / ystp[i]
        )
        # Subtracting new term from rhs term
        R[ivx] = (
            R[ivx] -
            bbottomx[j+1, 1] * etas[i+1, j+1] / ystpc[i+1] / ystp[i]
        )
    end
    return inz
end

end # module 