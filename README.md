# LaNMM_predictive_coding_paper

Code associated with CFC plots in Ruffini et al. 2025:
https://www.biorxiv.org/content/10.1101/2025.03.19.644090v3

LaNMM model reference:
https://www.sciencedirect.com/science/article/pii/S105381192300085X?via%3Dihub

## Repository layout

- `julia/`: current Julia implementation (engine, analysis, plots, sweep, run scripts).
- `python/`: legacy/original Python implementation (`lanmmv11` stack).
- `notebooks/`: exploratory notebooks and figure-building workflows.
- `julia_lanmm_sweep_*/`: generated sweep artifacts currently tracked in repo history.

## Quick start (2 minutes)

```bash
cd /path/to/LaNMM_predictive_coding_paper
julia --version
julia --threads auto julia/run_job_P1.jl
```

## Install Julia (macOS / Windows / Linux)

Recommended on all platforms: `juliaup` (official installer + version manager).

### macOS / Linux

```bash
curl -fsSL https://install.julialang.org | sh
julia --version
```

### Windows (PowerShell)

```powershell
winget install --id JuliaLang.Juliaup -e
julia --version
```

If `julia` is not found, restart your terminal/session.

## Run the Julia implementation

From the repository root:

```bash
julia --threads auto julia/run_job_P1.jl
julia --threads auto julia/run_job_P2.jl
julia --threads auto julia/run_job_PV.jl
```

- `julia/run_job_P1.jl`: multiscale driving on P1, constant on P2/PV.
- `julia/run_job_P2.jl`: multiscale driving on P2, constant on P1/PV.
- `julia/run_job_PV.jl`: multiscale driving on PV, constant on P1/P2.

## Julia architecture

- `julia/LaNMM_Types.jl`: typed config/result contracts.
- `julia/LaNMM_parameters.jl`: intrinsic constants and default driving model.
- `julia/LaNMM_Engine.jl`: ODE and single-simulation runtime.
- `julia/LaNMM_Analyzer.jl`: analysis metrics (SEC/EEC, power, PEIX, peaks).
- `julia/LaNMM_Plots.jl`: plotting only.
- `julia/LaNMM_Sweep.jl`: sweep orchestration, progress, saving.

Full Julia-specific walkthrough:

- `Readme_Julia_implemenation.md`

## Outputs per run

Each run writes a timestamped folder containing:

- `parameters.txt` (full reproducibility snapshot)
- diagnostic plots (`inputs.png`, `v.png`, `u_external.png`, `psd.png`)
- sweep plots (`couplings.png`, `power.png`, `peix.png`, `frequency_heatmaps.png`)
- serialized data (`sweep_results.jls`, `analysis_results.jls`)

## Notes

- First run is slower due to Julia package precompilation.
- For faster tests, reduce sweep ranges and `tmax` in `julia/run_job_*.jl`.
- Sweep output folders are currently versioned; very large `.jls` files can hit
  GitHub limits, so prefer archiving externally or using Git LFS for big runs.
