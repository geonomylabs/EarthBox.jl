module LoopTermination

function check_stop_global_loop_markers(
    iglobal::Int,
    nglobal::Int,
    convergence_criterion::Float64,
    tol_global::Float64
)::Bool
    check = false
    if iglobal == nglobal || convergence_criterion < tol_global
        check = true
    end
    return check
end

function check_stop_global_loop_nodes(
    nnode_yield::Int,
    iglobal::Int,
    nglobal::Int,
    global_yield_error::Float64,
    convergence_criterion::Float64,
    tolerance::Float64
)::Bool
    check = false
    if (nnode_yield == 0 || iglobal == nglobal ||
        global_yield_error < tolerance ||
        convergence_criterion < tolerance)
        check = true
    end
    return check
end

end # module 