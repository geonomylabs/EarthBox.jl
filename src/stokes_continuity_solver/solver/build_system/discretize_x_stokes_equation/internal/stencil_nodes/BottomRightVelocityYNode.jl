module BottomRightVelocityYNode

import EarthBox.Domain: basic_cell_along_upper_boundary
import EarthBox.Domain: basic_cell_along_lower_boundary
import ....BuildStructs: StokesBuildData
import ....BuildStructs: CellIndices

function coefficients_and_rhs_term(
    inz::Int, 
    inz_UR::Int, 
    inz_LR::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Int
    ivx = cell_indices.ivx
    i = cell_indices.i
    j = cell_indices.j
    ynum = build_data.grid.ynum
    xstpc = build_data.grid.xstpc
    ystp = build_data.grid.ystp
    etas = build_data.etas
    bbottomy = build_data.bc.bbottomy
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R
    hshift_to_vxR = build_data.grid.hshift_to_vxR

    if !basic_cell_along_lower_boundary(i, ynum)
        ivyLR = ivx + 1 + hshift_to_vxR
        if !basic_cell_along_upper_boundary(i)
            Lv[inz] = etas[i+1, j+1] / xstpc[j+1] / ystp[i]
            Li[inz] = ivx
            Lj[inz] = ivyLR
            Lii[inz] = i
            Ljj[inz] = j
            inz = inz + 1
        else
            # Basic node is at top boundary. Since coefficient was
            # already defined for vyLR we need to update the
            # previous value using the master array index.
            Lv[inz_LR] = (
                Lv[inz_LR] +
                etas[i+1, j+1] / xstpc[j+1] / ystp[i]
            )
        end
    else
        # Apply boundary condition since vyLR is a ghost node.
        # The boundary condition equation includes vyUR and,
        # therefore, a term is introduced that requires updating
        # the coefficient for vyUR. A constant is also introduced
        # that must be subtracted from the rhs term of the central node.
        Lv[inz_UR] = (
            Lv[inz_UR] +
            bbottomy[j+2, 2] * etas[i+1, j+1] / xstpc[j+1] / ystp[i]
        )
        R[ivx] = (
            R[ivx] -
            bbottomy[j+2, 1] * etas[i+1, j+1] / xstpc[j+1] / ystp[i]
        )
    end
    return inz
end

end # module 