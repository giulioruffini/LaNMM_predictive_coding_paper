using DSP
using Plots

"""
Plotting utilities for LaNMM Julia workflow.

Public API:
- `plot_sim_results`
- `plot_coupling_heatmaps`
- `plot_power_heatmaps_bottom_cbar`
- `plot_peix_heatmaps`
- `plot_frequency_heatmaps`
"""

function _matrix_from_sweep(dict_like, field)
    m1_vals = sort(unique(first(k) for k in keys(dict_like)))
    m2_vals = sort(unique(last(k) for k in keys(dict_like)))
    z = fill(NaN, length(m1_vals), length(m2_vals))
    for (i, m1) in enumerate(m1_vals), (j, m2) in enumerate(m2_vals)
        z[i, j] = getfield(dict_like[(m1, m2)], field)
    end
    return m1_vals, m2_vals, z
end

"""
Plot a single simulation's diagnostics (inputs, voltages, external u, PSD).
"""
function plot_sim_results(result::SimulationResult; save_dir::Union{Nothing,String}=nothing, xlims=nothing)
    n = minimum((
        length(result.t), length(result.e1_array), length(result.e2_array), length(result.pv_array),
        length(result.vP1), length(result.vP2), length(result.vPV),
        length(result.u3), length(result.u8), length(result.u14)
    ))
    if n < 4
        error("Not enough samples to plot (n=$n).")
    end
    t = result.t[1:n]
    e1 = result.e1_array[1:n]
    e2 = result.e2_array[1:n]
    pv = result.pv_array[1:n]
    vP1 = result.vP1[1:n]
    vP2 = result.vP2[1:n]
    vPV = result.vPV[1:n]
    u3 = result.u3[1:n]
    u8 = result.u8[1:n]
    u14 = result.u14[1:n]

    p_inputs = plot(t, e1, label="P1 input", color=:blue, linewidth=1.5,
        xlabel="Time (s)", ylabel="Firing Rate (Hz)", title="Inputs", grid=true)
    plot!(p_inputs, t, e2, label="P2 input", color=:red, linewidth=1.5)
    plot!(p_inputs, t, pv, label="PV input", color=:orange, linewidth=1.5)
    if xlims !== nothing
        xlims!(p_inputs, xlims...)
    end

    p_v = plot(t, vP1, label="vP1", color=:blue, linewidth=1.5,
        xlabel="Time (s)", ylabel="Membrane Potential (mV)",
        title="LaNMM Simulation Results (Membrane Potential)", grid=true)
    plot!(p_v, t, vP2, label="vP2", color=:red, linewidth=1.5)
    plot!(p_v, t, vPV, label="vPV", color=:orange, linewidth=1.5)
    if xlims !== nothing
        xlims!(p_v, xlims...)
    end

    p_uext = plot(t, u3, label="u3", color=:blue, linewidth=1.5,
        xlabel="Time (s)", ylabel="u (mV)", title="u Connected to External Drives", grid=true)
    plot!(p_uext, t, u8, label="u8", color=:red, linewidth=1.5)
    plot!(p_uext, t, u14, label="u14", color=:orange, linewidth=1.5)
    if xlims !== nothing
        xlims!(p_uext, xlims...)
    end

    p_psd = let
        combo = vP1 .+ vP2
        pgram = welch_pgram(combo, 2_048, 1_024; fs=result.fs)
        f = freq(pgram)
        s = max.(power(pgram), eps(Float64))
        plot(f, s, xlabel="Frequency (Hz)", ylabel="Power Spectral Density",
            title="Welch PSD of vP1+vP2", yscale=:log10, xlims=(0, 50), grid=true,
            linewidth=1.5, label="PSD")
    end

    if save_dir !== nothing
        mkpath(save_dir)
        savefig(p_inputs, joinpath(save_dir, "inputs.png"))
        savefig(p_v, joinpath(save_dir, "v.png"))
        savefig(p_uext, joinpath(save_dir, "u_external.png"))
        savefig(p_psd, joinpath(save_dir, "psd.png"))
    end

    display(plot(p_inputs, p_v, p_uext, p_psd, layout=(4, 1), size=(1200, 1200)))
    return nothing
end

function plot_coupling_heatmaps(couplings; title_str="Couplings", save_path=nothing)
    m1, m2, s2e = _matrix_from_sweep(couplings, :r_s2e)
    _, _, e2e = _matrix_from_sweep(couplings, :r_e2e)

    p1 = heatmap(m1, m2, s2e', color=:seismic, clims=(-0.7, 0.7),
        xlabel="P1 drive (Hz)", ylabel="P2 drive (Hz)", title="SEC Coupling", aspect_ratio=:equal)
    p2 = heatmap(m1, m2, e2e', color=:seismic, clims=(-0.7, 0.7),
        xlabel="P1 drive (Hz)", ylabel="P2 drive (Hz)", title="EEC Coupling", aspect_ratio=:equal)
    fig = plot(p1, p2, layout=(1, 2), size=(1200, 500), plot_title=title_str)
    if save_path !== nothing
        savefig(fig, save_path)
    end
    display(fig)
    return nothing
end

function plot_power_heatmaps_bottom_cbar(power_results; title_str="Band Power", save_path=nothing, log_scale=true)
    m1, m2, p1a = _matrix_from_sweep(power_results, :p1_alpha)
    _, _, p1g = _matrix_from_sweep(power_results, :p1_gamma)
    _, _, p2a = _matrix_from_sweep(power_results, :p2_alpha)
    _, _, p2g = _matrix_from_sweep(power_results, :p2_gamma)

    if log_scale
        floor_val = 1e-6
        p1a = log10.(max.(p1a, floor_val))
        p1g = log10.(max.(p1g, floor_val))
        p2a = log10.(max.(p2a, floor_val))
        p2g = log10.(max.(p2g, floor_val))
    end

    allv = vcat(vec(p1a), vec(p1g), vec(p2a), vec(p2g))
    vmin, vmax = minimum(allv), maximum(allv)

    kw = (xlabel="P1 drive (Hz)", ylabel="P2 drive (Hz)", color=:viridis, aspect_ratio=:equal, clims=(vmin, vmax))
    h1 = heatmap(m1, m2, p2a'; title="P2 Alpha Power", kw...)
    h2 = heatmap(m1, m2, p2g'; title="P2 Gamma Power", kw...)
    h3 = heatmap(m1, m2, p1a'; title="P1 Alpha Power", kw...)
    h4 = heatmap(m1, m2, p1g'; title="P1 Gamma Power", kw...)
    fig = plot(h1, h2, h3, h4, layout=(2, 2), size=(1200, 1000), plot_title=title_str)
    if save_path !== nothing
        savefig(fig, save_path)
    end
    display(fig)
    return nothing
end

function plot_peix_heatmaps(peix_results; title_str="PEIX", save_path=nothing)
    m1, m2, p1 = _matrix_from_sweep(peix_results, :peix_P1)
    _, _, p2 = _matrix_from_sweep(peix_results, :peix_P2)
    h1 = heatmap(m1, m2, p1', color=:jet, clims=(-1, 1), aspect_ratio=:equal,
        xlabel="P1 drive (Hz)", ylabel="P2 drive (Hz)", title="PEIX for P1")
    h2 = heatmap(m1, m2, p2', color=:jet, clims=(-1, 1), aspect_ratio=:equal,
        xlabel="P1 drive (Hz)", ylabel="P2 drive (Hz)", title="PEIX for P2")
    fig = plot(h1, h2, layout=(1, 2), size=(1200, 500), plot_title=title_str)
    if save_path !== nothing
        savefig(fig, save_path)
    end
    display(fig)
    return nothing
end

function plot_frequency_heatmaps(freq_data; save_path=nothing)
    h1 = heatmap(freq_data.m1_values, freq_data.m2_values, freq_data.alpha_peaks',
        color=:cool, clims=(6, 12), xlabel="P1 drive (Hz)", ylabel="P2 drive (Hz)",
        title="Peak Alpha Frequency in P1 (Hz)")
    h2 = heatmap(freq_data.m1_values, freq_data.m2_values, freq_data.gamma_peaks',
        color=:hot, clims=(30, 50), xlabel="P1 drive (Hz)", ylabel="P2 drive (Hz)",
        title="Peak Gamma Frequency in P2 (Hz)")
    fig = plot(h1, h2, layout=(1, 2), size=(1200, 500))
    if save_path !== nothing
        savefig(fig, save_path)
    end
    display(fig)
    return nothing
end
