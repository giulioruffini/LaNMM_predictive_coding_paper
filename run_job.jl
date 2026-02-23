include("LaNMM_Types.jl")
include("LaNMM_parameters.jl")
include("LaNMM_Engine.jl")
include("LaNMM_Analyzer.jl")
include("LaNMM_Plots.jl")
include("LaNMM_Sweep.jl")

# ==============================================================================
# 1. USER CONFIGURATION (PUBLIC CONTRACT)
# Set typed configs here, then execute with `run_sweep_job(...)`.
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
    mu_p1_values=collect(50.0:20.0:400.0),
    mu_p2_values=collect(50.0:20.0:400.0),
    tmax=10.0,   # change to 300.0 for production sweeps
    dt=0.001,
    discard=1.0,
    alpha_band=(8.0, 12.0),
    gamma_band=(30.0, 50.0),
    quiet_progress=false # true => fewer progress updates
)

# ==============================================================================
# 2. EXECUTE
# ==============================================================================
@time run_sweep_job(intrinsic_params, job_params)