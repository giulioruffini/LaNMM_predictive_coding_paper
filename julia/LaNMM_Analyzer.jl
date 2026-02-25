using DSP
using Statistics

"""
Public API (this file)
----------------------
Analysis:
- `compute_s2e_e2e`
- `analyze_sweep_couplings`
- `analyze_sweep_power`
- `analyze_peak_frequencies`
- `sweep_peix`
"""

"""
Bandpass filter with zero-phase forward-backward filtering.
"""
function bandpass_filter(signal, fs, lowcut, highcut; order=4)
    response = Bandpass(lowcut, highcut)
    filt = digitalfilter(response, Butterworth(order); fs=fs)
    return filtfilt(filt, signal)
end

"""
Analytic-signal envelope via Hilbert transform.
"""
get_envelope(signal) = abs.(hilbert(signal))

"""
Compute S2E and E2E coupling metrics from one simulation.

- S2E: corr(P1 alpha signal, P2 gamma envelope)
- E2E: corr(P1 alpha envelope, P2 gamma envelope)
"""
function compute_s2e_e2e(result; alpha_band=(8.0, 12.0), gamma_band=(30.0, 50.0), trim_seconds=0.5)
    t = result.t
    vP1 = result.vP1
    vP2 = result.vP2

    fs = 1.0 / (t[2] - t[1])
    p1_alpha = bandpass_filter(vP1, fs, alpha_band[1], alpha_band[2], order=4)
    p2_gamma = bandpass_filter(vP2, fs, gamma_band[1], gamma_band[2], order=4)

    env_p1 = get_envelope(p1_alpha)
    env_p2 = get_envelope(p2_gamma)

    trim_n = round(Int, trim_seconds * fs)
    if length(t) <= 2 * trim_n + 2
        return NaN, NaN
    end
    idx = (trim_n + 1):(length(t) - trim_n)

    r_s2e = cor(p1_alpha[idx], env_p2[idx])
    r_e2e = cor(env_p1[idx], env_p2[idx])
    return r_s2e, r_e2e
end

"""
Compute average band power with either Hilbert-envelope or Welch integration.
"""
function compute_band_power(signal, fs; freq_band=(8.0, 12.0), method=:bandpass_hilbert)
    if method == :bandpass_hilbert
        s = bandpass_filter(signal, fs, freq_band[1], freq_band[2], order=4)
        env = get_envelope(s)
        return mean(env .^ 2)
    elseif method == :welch
        pgram = welch_pgram(signal, 1024, 512; fs=fs)
        f = freq(pgram)
        p = power(pgram)
        idx = findall(x -> freq_band[1] <= x <= freq_band[2], f)
        return isempty(idx) ? 0.0 : sum(p[idx]) * (f[2] - f[1])
    end
    error("Unknown method: $method")
end

"""
Analyze S2E/E2E coupling over a sweep dictionary keyed by `(mu_p1, mu_p2)`.
"""
function analyze_sweep_couplings(sweep_results; alpha_band=(8.0, 12.0), gamma_band=(30.0, 50.0))
    # Snapshot keys once so threaded loop has deterministic indexable work.
    keys_vec = collect(keys(sweep_results))
    vals = Vector{NamedTuple{(:r_s2e, :r_e2e),Tuple{Float64,Float64}}}(undef, length(keys_vec))

    Threads.@threads for i in eachindex(keys_vec)
        key = keys_vec[i]
        result = sweep_results[key]
        r_s2e, r_e2e = compute_s2e_e2e(result; alpha_band=alpha_band, gamma_band=gamma_band)
        vals[i] = (r_s2e=r_s2e, r_e2e=r_e2e)
    end

    # Reconstruct dictionary only after all threads finished writing `vals`.
    return Dict(keys_vec[i] => vals[i] for i in eachindex(keys_vec))
end

"""
Analyze alpha/gamma power over the sweep for P1 and P2.
"""
function analyze_sweep_power(sweep_results; alpha_band=(8.0, 12.0), gamma_band=(30.0, 50.0), method=:bandpass_hilbert)
    keys_vec = collect(keys(sweep_results))
    vals = Vector{NamedTuple{(:p1_alpha,:p1_gamma,:p2_alpha,:p2_gamma),Tuple{Float64,Float64,Float64,Float64}}}(undef, length(keys_vec))

    Threads.@threads for i in eachindex(keys_vec)
        key = keys_vec[i]
        result = sweep_results[key]
        fs = result.fs
        vals[i] = (
            p1_alpha=compute_band_power(result.vP1, fs; freq_band=alpha_band, method=method),
            p1_gamma=compute_band_power(result.vP1, fs; freq_band=gamma_band, method=method),
            p2_alpha=compute_band_power(result.vP2, fs; freq_band=alpha_band, method=method),
            p2_gamma=compute_band_power(result.vP2, fs; freq_band=gamma_band, method=method)
        )
    end

    return Dict(keys_vec[i] => vals[i] for i in eachindex(keys_vec))
end

"""
Peak frequency in `band` using Welch spectrum.
"""
function compute_peak_frequency(signal, fs, band)
    pgram = welch_pgram(signal, 1024, 512; fs=fs)
    freqs = freq(pgram)
    p = power(pgram)
    idx = findall(f -> band[1] <= f <= band[2], freqs)
    if isempty(idx)
        return 0.0
    end
    local_idx = argmax(p[idx])
    return freqs[idx[local_idx]]
end

"""
Analyze alpha (P1) and gamma (P2) peak-frequency maps over the sweep.
"""
function analyze_peak_frequencies(sweep_results; fs=1000.0, alpha_band=(8.0, 12.0), gamma_band=(30.0, 80.0))
    m1_vals = sort(unique(first(k) for k in keys(sweep_results)))
    m2_vals = sort(unique(last(k) for k in keys(sweep_results)))

    alpha_peaks = zeros(length(m1_vals), length(m2_vals))
    gamma_peaks = zeros(length(m1_vals), length(m2_vals))

    total = length(m1_vals) * length(m2_vals)
    Threads.@threads for idx in 1:total
        i = ((idx - 1) % length(m1_vals)) + 1
        j = ((idx - 1) ÷ length(m1_vals)) + 1
        m1 = m1_vals[i]
        m2 = m2_vals[j]
        r = sweep_results[(m1, m2)]
        alpha_peaks[i, j] = compute_peak_frequency(r.vP1, fs, alpha_band)
        gamma_peaks[i, j] = compute_peak_frequency(r.vP2, fs, gamma_band)
    end

    return (alpha_peaks=alpha_peaks, gamma_peaks=gamma_peaks, m1_values=m1_vals, m2_values=m2_vals)
end

"""
PEIX nonlinearity helper.
"""
compute_peix(x) = -sign(x) * (4.0 * exp(-x) / (1.0 + exp(-x))^2 - 1.0)

"""
Sweep PEIX values for P1 and P2 from sweep simulation outputs.
"""
function sweep_peix(sweep_results; r_slope=0.56, v0_default=6.0, v0_p2=1.0)
    keys_vec = collect(keys(sweep_results))
    vals = Vector{NamedTuple{(:peix_P1,:peix_P2),Tuple{Float64,Float64}}}(undef, length(keys_vec))

    Threads.@threads for i in eachindex(keys_vec)
        key = keys_vec[i]
        result = sweep_results[key]
        x1 = r_slope * (mean(result.vP1) - v0_default)
        x2 = r_slope * (mean(result.vP2) - v0_p2)
        vals[i] = (peix_P1=compute_peix(x1), peix_P2=compute_peix(x2))
    end

    return Dict(keys_vec[i] => vals[i] for i in eachindex(keys_vec))
end