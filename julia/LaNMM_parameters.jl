#=
Column-level intrinsic LaNMM parameters and condition modifiers.

This file centralizes constants that define the local cortical microcircuit:
- synaptic gains/time constants by synapse type
- synapse type assignment for the 14 synapses
- baseline connectivity constants (C values)
- condition-dependent connectivity modifiers
=#
if !isdefined(@__MODULE__, :DrivingConfig)
    include("LaNMM_Types.jl")
end

"""
Public API (this file)
----------------------
- `get_column_parameters(; condition="healthy")`
    Return intrinsic column constants with condition modifiers applied.
- `build_synapse_parameter_arrays(column_params)`
    Convert column constants into ODE-ready `(a_vals, A_vals, C_vals)`.
- `get_default_driving_model()`
    Return canonical driving defaults used by the current Julia engine.
"""

const LANMM_PSY_FACTOR = 1.3846153846
const LANMM_SUPPORTED_CONDITIONS = ("healthy", "mci", "ad", "psychedelics", "mci+psy", "ad+psy")

const LANMM_BASE_CONNECTIVITY = Float64[
    108.0,  # C1
    33.7,   # C2
    1.0,    # C3  (external e1)
    135.0,  # C4
    33.75,  # C5
    70.0,   # C6
    550.0,  # C7
    1.0,    # C8  (external e2)
    200.0,  # C9
    100.0,  # C10
    80.0,   # C11
    200.0,  # C12
    30.0,   # C13
    1.0     # C14 (external pv)
]

const LANMM_SYNAPSE_TYPES = (
    :AMPA, :GABA_slow, :AMPA, :AMPA, :AMPA, :AMPA, :GABA_fast,
    :AMPA, :AMPA, :GABA_fast, :AMPA, :AMPA, :AMPA, :AMPA
)

const LANMM_DEFAULT_DRIVING_MODEL = DrivingConfig(
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

"""
Return column-level parameters after applying condition-specific modifiers.
"""
function get_column_parameters(; condition::String="healthy")
    cond = lowercase(condition)
    if !(cond in LANMM_SUPPORTED_CONDITIONS)
        error("Unsupported condition '$condition'. Supported: $(collect(LANMM_SUPPORTED_CONDITIONS))")
    end

    params = (
        A_AMPA=3.25,
        a_AMPA=100.0,
        A_GABA_slow=-22.0,
        a_GABA_slow=50.0,
        A_GABA_fast=-30.0,
        a_GABA_fast=220.0,
        C_vals=copy(LANMM_BASE_CONNECTIVITY),
        syn_types=LANMM_SYNAPSE_TYPES,
        # Sigmoid defaults are kept here so the same source can be reused
        # by scripts that want a canonical intrinsic setup.
        v0_default=6.0,
        v0_p2=1.0,
        fmax=5.0,
        r_slope=0.56
    )

    if cond == "mci" || cond == "mci+psy"
        params.C_vals[7] = 300.0
    elseif cond == "ad" || cond == "ad+psy"
        params.C_vals[7] = 140.0
    end

    if occursin("psy", cond) || cond == "psychedelics"
        for syn in (1, 3, 11)
            params.C_vals[syn] *= LANMM_PSY_FACTOR
        end
    end

    return params
end

"""
Return canonical driving defaults used by the Julia engine and sweep metadata.
"""
get_default_driving_model() = LANMM_DEFAULT_DRIVING_MODEL

"""
Build fast per-synapse vectors `(a_vals, A_vals, C_vals)` for the ODE.
"""
function build_synapse_parameter_arrays(column_params)
    A_vals = zeros(Float64, 14)
    a_vals = zeros(Float64, 14)
    C_vals = copy(column_params.C_vals)

    for i in 1:14
        st = column_params.syn_types[i]
        if st == :AMPA
            A_vals[i] = column_params.A_AMPA
            a_vals[i] = column_params.a_AMPA
        elseif st == :GABA_slow
            A_vals[i] = column_params.A_GABA_slow
            a_vals[i] = column_params.a_GABA_slow
        else
            A_vals[i] = column_params.A_GABA_fast
            a_vals[i] = column_params.a_GABA_fast
        end
    end

    return a_vals, A_vals, C_vals
end
