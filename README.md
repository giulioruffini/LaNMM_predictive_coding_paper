# LaNMM_predictive_coding_paper
Code associated with CFC plots in Ruffini et al 2025 https://www.biorxiv.org/content/10.1101/2025.03.19.644090v3
(Cross-Frequency Coupling as a Neural Substrate for Prediction Error Evaluation: A Laminar Neural Mass Modeling Approach, 2025)

The LaNMM model is as the one described in https://www.sciencedirect.com/science/article/pii/S105381192300085X?via%3Dihub 
(A physical neural mass model framework for the analysis of oscillatory generators from laminar electrophysiological recordings, Sanchez-Todo et al., 2023)

Notes: 
The LaNMM is also used  https://www.biorxiv.org/content/10.1101/2024.12.15.628565v4
(Restoring Oscillatory Dynamics in Alzheimer’s Disease: A Laminar Whole-Brain Model of Serotonergic Psychedelic Effects, 2024)
In the code, we chose to modify the connectivity constants ($C$) rather than the global synaptic gain ($A$). 
This was done to achieve specificity in a simple manner (a hack, so to speak): modifying $A$ would affect 
all excitatory synapses globally, forcing us to rewrite the code for more analytic control of the variable at all synapses, 
whereas modifying $C$ allowed us to selectively target the specific glutamatergic inputs to the P1 population as intended. 
We have added comments to the code repository to explicitly state this equivalence: scaling $C$ 
for specific synapses is our implementation method for increasing the effective synaptic gain for those specific pathways.



The LaNMM is also used in https://www.biorxiv.org/content/10.1101/2025.03.26.645407v1
(Fast Interneuron Dysfunction in Laminar Neural Mass Model Reproduces Alzheimer’s Oscillatory Biomarkers, 2025)

(and probably more to follow)

## Julia workflow (parity with Python plots)

The Julia code is split into:

- `LaNMM_Types.jl`: typed configs and simulation result contract.
- `LaNMM_parameters.jl`: intrinsic column constants + driving defaults.
- `LaNMM_Engine.jl`: core ODE model + simulation (`run_unified_simulation`).
- `LaNMM_Analyzer.jl`: analysis functions only (no plotting side effects).
- `LaNMM_Plots.jl`: plotting functions only.
- `LaNMM_Sweep.jl`: sweep orchestration/progress/output helpers.
- `run_job.jl`: clean entry script to define typed parameters and execute.

Default Julia driving setup:
- `e1`: multiscale
- `e2`: constant
- `pv`: constant

Run:

`julia --threads auto run_job.jl`

Detailed Julia architecture and API documentation:

- `Readme_Julia_implemenation.md`

Each run creates a timestamped output directory containing:

- `inputs.png`, `v.png`, `u_external.png`, `psd.png` (single-run diagnostics)
- `couplings.png` (SEC/EEC heatmaps)
- `power.png` (P1/P2 alpha/gamma band-power heatmaps)
- `peix.png` (PEIX heatmaps)
- `frequency_heatmaps.png` (alpha/gamma peak-frequency heatmaps)
- `sweep_results.jls`, `analysis_results.jls`, and `parameters.txt`

### Python-to-Julia function mapping

- Python `run_unified_simulation` -> Julia `run_unified_simulation` (`LaNMM_Engine.jl`)
- Python `plot_sim_results` -> Julia `plot_sim_results` (`LaNMM_Analyzer.jl`)
- Python `analyze_sweep_couplings` / `plot_coupling_heatmaps` -> same names in Julia
- Python `analyze_sweep_power` / `plot_power_heatmaps_bottom_cbar` -> same names in Julia
- Python `sweep_peix` / `plot_peix_heatmaps` -> same names in Julia
- Python `analyze_peak_frequencies` / `plot_frequency_heatmaps` -> same names in Julia
