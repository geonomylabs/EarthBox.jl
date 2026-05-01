module RunitLoop

using EarthBox

function run_local_models()::Nothing
    base_path = "/home/$(get_username())/apps/julia/EarthBox/models/extension_to_sfs"
    models = Dict(
        "$(base_path)/extension_strong_zones" => ["case1"],
        "$(base_path)/extension_weak_fault" => ["case2"],
    )
    local_model_loop(
        models       = models,
        run_model    = true,
        plot_markers = false,
        plot_scalars = false,
        istart       = 1,
        iend         = 100
    )
    return nothing
end

function main()
    run_local_models()
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    RunitLoop.main()
end
