module BackBottomVelocityYNode

import EarthBox.Domain: basic_cell_along_upper_boundary
import EarthBox.Domain: basic_cell_along_lower_boundary
import ....BuildStructs: StokesBuildData
import ....BuildStructs: CellIndices

function coefficients_and_rhs_term(
    inz::Int,
    inz_BU::Int,
    inz_BD::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Int
    i = cell_indices.i
    j = cell_indices.j
    k = cell_indices.k
    ivz = cell_indices.ivz
    ynum = build_data.grid.ynum
    ystp = build_data.grid.ystp
    zstpc = build_data.grid.zstpc
    etas_yz = build_data.etas_yz
    bbottomy = build_data.bc.bbottomy
    dshift_to_vxF = build_data.grid.dshift_to_vxF
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    if !basic_cell_along_lower_boundary(i, ynum)
        ivyBD = ivz - 1 + dshift_to_vxF
        if !basic_cell_along_upper_boundary(i)
            Lv[inz] = etas_yz[i+1, j, k+1] / ystp[i] / zstpc[k+1]
            Li[inz] = ivz
            Lj[inz] = ivyBD
            Lii[inz] = i
            Ljj[inz] = j
            inz = inz + 1
        else
            # Basic cell is at top boundary; vyBD coefficient was already
            # defined by BackUpperVelocityYNode at inz_BD. Update it.
            Lv[inz_BD] = (
                Lv[inz_BD] +
                etas_yz[i+1, j, k+1] / ystp[i] / zstpc[k+1]
            )
        end
    else
        # vyBD at bottom boundary. Substitute via bbottomy.
        Lv[inz_BU] = (
            Lv[inz_BU] +
            bbottomy[j+1, k+2, 2] * etas_yz[i+1, j, k+1] / ystp[i] / zstpc[k+1]
        )
        R[ivz] = (
            R[ivz] -
            bbottomy[j+1, k+2, 1] * etas_yz[i+1, j, k+1] / ystp[i] / zstpc[k+1]
        )
    end
    return inz
end

end # module
