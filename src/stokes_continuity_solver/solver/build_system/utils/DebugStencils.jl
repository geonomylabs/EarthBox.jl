module DebugStencils

import ..BuildStructs: CellIndices
import ..BuildStructs: StokesBuildData

const TARGET_INZ = 30474

function check_inz_value_continuity_internal(inz::Int64, loc_id::Int64)::Nothing
    if inz == TARGET_INZ
        println("Continuity Internal $loc_id: inz = $inz")
    end
    return nothing
end

function check_inz_value_continuity_boundary(inz::Int64, loc_id::Int64)::Nothing
    if inz == TARGET_INZ
        println("Continuity Boundary $loc_id: inz = $inz")
    end
    return nothing
end

function check_x_stokes_internal(
    inz::Int64, 
    loc_id::Int64, 
    cell_indices::CellIndices, 
    build_data::StokesBuildData
)::Nothing
    check_inz_value_x_stokes_internal(inz, loc_id)
    check_large_matrix_x_stokes(inz, loc_id, cell_indices, build_data)
    return nothing
end

function check_inz_value_x_stokes_internal(inz::Int64, loc_id::Int64, )::Nothing
    if inz == TARGET_INZ
        println("X-Stokes Internal $loc_id: inz = $inz")
    end
    return nothing
end

function check_large_matrix_x_stokes(
    inz::Int64, 
    loc_id::Int64, 
    cell_indices::CellIndices, 
    build_data::StokesBuildData
)::Nothing
    if inz == TARGET_INZ
        Lv = build_data.system_vectors.Lv[inz]
        Lii = build_data.system_vectors.Lii[inz]
        Ljj = build_data.system_vectors.Ljj[inz]
        i = cell_indices.i
        j = cell_indices.j
        ivx = cell_indices.ivx
        ivy = cell_indices.ivy
        ipr = cell_indices.ipr
        println("Large Matrix X-Stokes $loc_id: inz = $inz, Lv = $Lv, Lii = $Lii, Ljj = $Ljj, i = $i, j = $j, ivx = $ivx, ivy = $ivy, ipr = $ipr")
    end
    return nothing
end

function check_inz_value_x_stokes_boundary(inz::Int64, loc_id::Int64)::Nothing
    if inz == TARGET_INZ
        println("X-Stokes Boundary $loc_id: inz = $inz")
    end
    return nothing
end

function check_inz_value_y_stokes_internal(inz::Int64, loc_id::Int64)::Nothing
    if inz == TARGET_INZ
        println("Y-Stokes Internal $loc_id: inz = $inz")
    end
    return nothing
end

function check_inz_value_y_stokes_boundary(inz::Int64, loc_id::Int64)::Nothing
    if inz == TARGET_INZ
        println("Y-Stokes Boundary $loc_id: inz = $inz")
    end
    return nothing
end

function check_large_matrix_in_stencil(
    inz::Int64,
    name::String,
    i::Int64, 
    j::Int64,
    Lv::Vector{Float64},
    etan::Float64,
    xstp::Float64,
    xstpc::Float64
)::Nothing
    if inz == TARGET_INZ
        println(" $name inz = $inz, Lv[$inz] = $(Lv[inz]), i = $i, j = $j, etan = $etan, xstp = $xstp, xstpc = $xstpc")
    end
    return nothing
end

end