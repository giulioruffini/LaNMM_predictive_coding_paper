using DifferentialEquations
using DataInterpolations
using Statistics

"""
Container for all parameters required by the ODE function.

The interpolation fields are parametric so the struct remains type-stable even if
the interpolation implementation changes.
"""
Base.@kwdef struct LaNMMParams{T1,T2,T3}
    v0_default::Float64
    v0_p2::Float64
    fmax::Float64
    r_slope::Float64
    a_vals::Vector{Float64}
    A_vals::Vector{Float64}
    C_vals::Vector{Float64}
    e1_interp::T1
    e2_interp::T2
    pv_interp::T3
end

"""
Build intrinsic synaptic arrays for a given condition.

Returns `(a_vals, A_vals, C_vals)` for synapses 1..14.
"""
function build_intrinsic_params(; condition::String="healthy")
    A_AMPA = 3.25
    a_AMPA = 100.0
    A_GABA_slow = -22.0
    a_GABA_slow = 50.0
    A_GABA_fast = -30.0
    a_GABA_fast = 220.0

    C_dict = Dict(
        1 => 108.0, 2 => 33.7, 3 => 1.0, 4 => 135.0, 5 => 33.75, 6 => 70.0,
        7 => 550.0, 8 => 1.0, 9 => 200.0, 10 => 100.0, 11 => 80.0, 12 => 200.0,
        13 => 30.0, 14 => 1.0
    )
    syn_types = Dict(
        1 => "AMPA", 2 => "GABA_slow", 3 => "AMPA", 4 => "AMPA", 5 => "AMPA",
        6 => "AMPA", 7 => "GABA_fast", 8 => "AMPA", 9 => "AMPA", 10 => "GABA_fast",
        11 => "AMPA", 12 => "AMPA", 13 => "AMPA", 14 => "AMPA"
    )

    cond = lowercase(condition)
    if cond == "mci" || cond == "mci+psy"
        C_dict[7] = 300.0
    elseif cond == "ad" || cond == "ad+psy"
        C_dict[7] = 140.0
    end

    if occursin("psy", cond) || cond == "psychedelics"
        for syn in (1, 3, 11)
            C_dict[syn] *= 1.3846153846
        end
    end

    A_vals = zeros(Float64, 14)
    a_vals = zeros(Float64, 14)
    C_vals = zeros(Float64, 14)

    for i in 1:14
        C_vals[i] = C_dict[i]
        if syn_types[i] == "AMPA"
            A_vals[i] = A_AMPA
            a_vals[i] = a_AMPA
        elseif syn_types[i] == "GABA_slow"
            A_vals[i] = A_GABA_slow
            a_vals[i] = a_GABA_slow
        else
            A_vals[i] = A_GABA_fast
            a_vals[i] = a_GABA_fast
        end
    end

    return a_vals, A_vals, C_vals
end

"""
Generate a simple amplitude-modulated input signal centered around `baseline`.
"""
function generate_am_signal(t_array::AbstractVector, carrier_freq, carrier_amplitude, baseline)
    envelope = carrier_amplitude .* (1.0 .+ 0.5 .* sin.(2π .* 10.0 .* t_array))
    carrier = cos.(2π .* carrier_freq .* t_array)
    s = envelope .* carrier
    return (s .- mean(s)) .+ baseline
end

@inline sigmoid_pop(v, v0, p::LaNMMParams) = p.fmax / (1.0 + exp(p.r_slope * (v0 - v)))

"""
LaNMM right-hand side with 14 second-order synapses (28 states).
"""
function lanmm_ode!(du, u, p::LaNMMParams, t)
    phi_e1 = max(0.0, p.e1_interp(t))
    phi_e2 = max(0.0, p.e2_interp(t))
    phi_pv = max(0.0, p.pv_interp(t))

    @inbounds begin
        vP1 = u[1] + u[3] + u[5] + u[21]
        vSS = u[7]
        vSST = u[9]
        vP2 = u[11] + u[13] + u[15] + u[23]
        vPV = u[17] + u[19] + u[25] + u[27]

        sP1 = sigmoid_pop(vP1, p.v0_default, p)
        sSS = sigmoid_pop(vSS, p.v0_default, p)
        sSST = sigmoid_pop(vSST, p.v0_default, p)
        sP2 = sigmoid_pop(vP2, p.v0_p2, p)
        sPV = sigmoid_pop(vPV, p.v0_default, p)

        presyn = (
            sSS, sSST, phi_e1, sP1, sP1, sP2, sPV,
            phi_e2, sP2, sPV, sP2, sP1, sP1, phi_pv
        )

        for s in 1:14
            idx_u = 2s - 1
            idx_z = 2s
            du[idx_u] = u[idx_z]
            du[idx_z] = p.a_vals[s] * p.A_vals[s] * (p.C_vals[s] * presyn[s]) -
                        2.0 * p.a_vals[s] * u[idx_z] -
                        (p.a_vals[s]^2) * u[idx_u]
        end
    end
    return nothing
end

"""
Run the unified Julia simulation and return a Python-like result payload.

Returned fields match the analyzer expectations:
`t`, `vP1`, `vP2`, `vPV`, `e1_array`, `e2_array`, `pv_array`, `fs`, `u1..u14`.
"""
function run_unified_simulation(;
    condition::String="healthy",
    mu_e1::Float64=270.0,
    mu_e2::Float64=90.0,
    mu_pv::Float64=0.0,
    v0_default::Float64=6.0,
    v0_p2::Float64=1.0,
    fmax::Float64=5.0,
    r_slope::Float64=0.56,
    tmax::Float64=4.0,
    dt::Float64=0.001,
    discard::Float64=1.0
)
    t_array = collect(0:dt:tmax)
    a_vals, A_vals, C_vals = build_intrinsic_params(condition=condition)

    e1_array = max.(0.0, generate_am_signal(t_array, 10.0, 400.0, mu_e1))
    e2_array = max.(0.0, generate_am_signal(t_array, 40.0, 400.0, mu_e2))
    pv_array = fill(mu_pv, length(t_array))

    p = LaNMMParams(
        v0_default=v0_default,
        v0_p2=v0_p2,
        fmax=fmax,
        r_slope=r_slope,
        a_vals=a_vals,
        A_vals=A_vals,
        C_vals=C_vals,
        e1_interp=LinearInterpolation(e1_array, t_array),
        e2_interp=LinearInterpolation(e2_array, t_array),
        pv_interp=LinearInterpolation(pv_array, t_array)
    )

    u0 = zeros(Float64, 28)
    prob = ODEProblem(lanmm_ode!, u0, (0.0, tmax), p)
    sol = solve(prob, Tsit5(), saveat=dt, reltol=1e-8, abstol=1e-8, maxiters=10^8)
    if !SciMLBase.successful_retcode(sol)
        # Fallback to a stiff-capable solver for hard parameter combinations.
        sol = solve(prob, Rosenbrock23(), saveat=dt, reltol=1e-8, abstol=1e-8, maxiters=10^8)
    end
    if !SciMLBase.successful_retcode(sol)
        error("ODE solve failed with retcode=$(sol.retcode). Try reducing dt/tmax or switching condition/drive settings.")
    end

    idx_start = searchsortedfirst(sol.t, discard)
    if idx_start > length(sol.t)
        error("Discard time ($discard s) is beyond the simulated interval ending at $(sol.t[end]) s.")
    end
    t_out = sol.t[idx_start:end]
    n = length(t_out)

    u_mat = Matrix{Float64}(undef, n, 28)
    for (k, uk) in enumerate(sol.u[idx_start:end])
        @inbounds u_mat[k, :] = uk
    end

    # State mapping follows the Python model convention:
    # u(s) corresponds to first state of synapse s (odd positions in the 28-state vector).
    u_syn(s) = @view u_mat[:, 2s - 1]
    vP1 = u_syn(1) .+ u_syn(2) .+ u_syn(3) .+ u_syn(11)
    vP2 = u_syn(6) .+ u_syn(7) .+ u_syn(8) .+ u_syn(12)
    vPV = u_syn(9) .+ u_syn(10) .+ u_syn(13) .+ u_syn(14)

    # Re-evaluate drives exactly on simulated output times so plotting arrays
    # stay aligned even if the integrator terminates early.
    e1_out = [p.e1_interp(tt) for tt in t_out]
    e2_out = [p.e2_interp(tt) for tt in t_out]
    pv_out = [p.pv_interp(tt) for tt in t_out]

    return (
        t=t_out,
        vP1=collect(vP1),
        vP2=collect(vP2),
        vPV=collect(vPV),
        e1_array=e1_out,
        e2_array=e2_out,
        pv_array=pv_out,
        fs=1.0 / dt,
        u1=collect(u_syn(1)),
        u2=collect(u_syn(2)),
        u3=collect(u_syn(3)),
        u4=collect(u_syn(4)),
        u5=collect(u_syn(5)),
        u6=collect(u_syn(6)),
        u7=collect(u_syn(7)),
        u8=collect(u_syn(8)),
        u9=collect(u_syn(9)),
        u10=collect(u_syn(10)),
        u11=collect(u_syn(11)),
        u12=collect(u_syn(12)),
        u13=collect(u_syn(13)),
        u14=collect(u_syn(14))
    )
end