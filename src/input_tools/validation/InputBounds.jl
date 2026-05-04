module InputBounds

# DoS-prevention bounds on model.yml fields. Numbers are chosen well above any
# realistic scientific configuration so legitimate inputs pass unchanged. They
# exist to fail fast on configs that would otherwise OOM the host or run loops
# for years.

const INT_BOUNDS = Dict{String, Tuple{Int, Int}}(
    "xnum"            => (1, 1_000_000),
    "ynum"            => (1, 1_000_000),
    "ntimestep_max"   => (0, 10_000_000),
    "nglobal"         => (0, 10_000),
    "nsand_layers"    => (0, 1_000),
    "ntime_increase_1" => (0, 10_000_000),
    "ntime_increase_2" => (0, 10_000_000),
    "ntime_increase_3" => (0, 10_000_000),
    "ntime_increase_4" => (0, 10_000_000),
)

# Float fields where the value drives allocation or counts. Treated as Real so
# they accept Int or Float input.
const REAL_BOUNDS = Dict{String, Tuple{Real, Real}}(
    "nmarkers_cell_x" => (1, 100),
    "nmarkers_cell_y" => (1, 100),
)

# Float fields that must be finite. Adding a name here means NaN/Inf reject.
const FLOAT_FINITE = Set{String}([
    "xsize", "ysize",
    "viscosity_min", "viscosity_max",
    "yield_stress_min", "yield_stress_max",
    "timestep_viscoelastic", "timestep_out",
    "displ_limit",
    "max_temp_change",
    "tolerance_picard",
    "memory_relax_perc",
    "pressure_bc",
    "temperature_top", "temperature_bottom",
    "gravity_x", "gravity_y",
    "marker_cell_displ_max",
    "subgrid_diff_coef_stress", "subgrid_diff_coef_temp",
    "cell_displ_factor", "time_increase_factor",
])

# Float fields that must additionally be > 0.
const FLOAT_POSITIVE = Set{String}([
    "xsize", "ysize",
    "viscosity_min", "viscosity_max",
    "timestep_viscoelastic", "timestep_out",
    "max_temp_change",
])

# Total markers cap (xnum * ynum * nmcx * nmcy). Per-field caps alone permit
# 10^16 in the worst case, so a derived check is required.
const TOTAL_MARKERS_CAP = 10_000_000_000  # 10^10

"""
    validate_value(name::String, value)

Apply bounds and finiteness checks for a known input parameter. Unknown names
pass silently — Reader.jl already rejects unknown names earlier in its own
registry check, so duplicating that logic here would only produce noise.

Throws `ArgumentError` with a message identifying the parameter and the
violated rule.
"""
function validate_value(name::String, value)::Nothing
    if haskey(INT_BOUNDS, name)
        check_int_bounds(name, value)
    end
    if haskey(REAL_BOUNDS, name)
        check_real_bounds(name, value)
    end
    if name in FLOAT_FINITE
        check_finite(name, value)
    end
    if name in FLOAT_POSITIVE
        check_positive(name, value)
    end
    return nothing
end

function check_int_bounds(name::String, value)::Nothing
    if !(value isa Integer)
        throw(ArgumentError(
            "Input parameter '$name' must be an integer; got $(typeof(value))."
        ))
    end
    lo, hi = INT_BOUNDS[name]
    if value < lo || value > hi
        throw(ArgumentError(
            "Input parameter '$name' = $value is outside the allowed range [$lo, $hi]."
        ))
    end
    return nothing
end

function check_real_bounds(name::String, value)::Nothing
    if !(value isa Real)
        throw(ArgumentError(
            "Input parameter '$name' must be a real number; got $(typeof(value))."
        ))
    end
    if value isa AbstractFloat && !isfinite(value)
        throw(ArgumentError(
            "Input parameter '$name' must be finite; got $value."
        ))
    end
    lo, hi = REAL_BOUNDS[name]
    if value < lo || value > hi
        throw(ArgumentError(
            "Input parameter '$name' = $value is outside the allowed range [$lo, $hi]."
        ))
    end
    return nothing
end

function check_finite(name::String, value)::Nothing
    if value isa AbstractFloat && !isfinite(value)
        throw(ArgumentError(
            "Input parameter '$name' must be finite; got $value."
        ))
    end
    return nothing
end

function check_positive(name::String, value)::Nothing
    if value isa Real && value <= 0
        throw(ArgumentError(
            "Input parameter '$name' must be > 0; got $value."
        ))
    end
    return nothing
end

"""
    validate_total_markers(xnum, ynum, nmarkers_cell_x, nmarkers_cell_y)

Guard against grid configurations that, while individually within per-field
bounds, multiply out to an unreasonable total marker count. Throws
`ArgumentError` if `xnum * ynum * nmarkers_cell_x * nmarkers_cell_y` exceeds
`TOTAL_MARKERS_CAP` (10^10).
"""
function validate_total_markers(xnum, ynum, nmarkers_cell_x, nmarkers_cell_y)::Nothing
    total = Float64(xnum) * Float64(ynum) *
            Float64(nmarkers_cell_x) * Float64(nmarkers_cell_y)
    if total > TOTAL_MARKERS_CAP
        throw(ArgumentError(
            "Total marker count $total = xnum × ynum × nmarkers_cell_x × " *
            "nmarkers_cell_y exceeds the allowed cap $TOTAL_MARKERS_CAP."
        ))
    end
    return nothing
end

end # module InputBounds
