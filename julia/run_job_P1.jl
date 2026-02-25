include("LaNMM_Types.jl")
include("LaNMM_parameters.jl")
include("LaNMM_Engine.jl")
include("LaNMM_Analyzer.jl")
include("LaNMM_Plots.jl")
include("LaNMM_Sweep.jl")

# ==============================================================================
# 1. USER CONFIGURATION (START HERE)
#
# New Julia users: this is the only section you usually need to edit.
# - `intrinsic_params` controls biological/intrinsic assumptions.
# - `job_params` controls sweep size and simulation timing.
# - `driving_params` chooses which population gets multiscale input.
# ==============================================================================
const intrinsic_params = IntrinsicConfig(
    condition="healthy",
    v0_default=6.0,
    v0_p2=1.0,
    fmax=5.0,
    r_slope=0.56
)

const job_params = JobConfig(
    job_title="julia_lanmm_sweep",
    mu_p1_values=collect(50.0:5.0:400.0),
    mu_p2_values=collect(50.0:5.0:400.0),
    tmax=603.0,   # change to 300.0 for production sweeps
    dt=0.001,
    discard=3.0,
    alpha_band=(8.0, 12.0),
    gamma_band=(30.0, 50.0),
    quiet_progress=false # true => fewer progress updates
)

# P1 multiscale, P2/PV constant
const driving_params = DrivingConfig(
    e1=MultiscaleDrivingConfig(
        mode=:multiscale,
        seed_offset=0,
        slow_std=400.0,
        slow_alpha=0.99,
        fast_std=5.0,
        fast_cutoff=100.0,
        floor_val=1e-4
    ),
    e2=ConstantDrivingConfig(mode=:constant, baseline=0.0),
    pv=ConstantDrivingConfig(mode=:constant, baseline=0.0),
    seed=42,
    nonnegative_clipping=true
)

# ==============================================================================
# 2. EXECUTE
# ==============================================================================
# This launches one full sweep, saves plots/data, and prints progress + ETA.
@time run_sweep_job(intrinsic_params, job_params; driving=driving_params)