# Readme Julia Implementation

This document describes the Julia implementation architecture, file responsibilities,
data contracts, and execution flow for the LaNMM pipeline.

All file paths below are relative to the `julia/` folder.

## Architecture Overview

The Julia code is organized by concern:

- `LaNMM_Types.jl`  
  Typed data contracts for configs and outputs.
- `LaNMM_parameters.jl`  
  Intrinsic column parameters and driving defaults.
- `LaNMM_Engine.jl`  
  ODE core and numerical simulation.
- `LaNMM_Analyzer.jl`  
  Analysis functions only (coupling, power, PEIX, peaks).
- `LaNMM_Plots.jl`  
  Plotting functions only.
- `LaNMM_Sweep.jl`  
  Sweep orchestration, progress/ETA, metadata snapshot, persistence.
- `run_job_P1.jl`, `run_job_P2.jl`, `run_job_PV.jl`  
  User-facing entry scripts: define configs, execute one run.

This separation keeps model equations, analysis logic, plotting, and job orchestration
independent and easier to reason about.

## Typed Contracts (Phase 2)

`LaNMM_Types.jl` introduces typed interfaces:

- `IntrinsicConfig`  
  `condition`, sigmoid parameters (`v0_default`, `v0_p2`, `fmax`, `r_slope`).
- `JobConfig`  
  sweep ranges, simulation timing (`tmax`, `dt`, `discard`), analysis bands, and progress verbosity.
- `DrivingConfig` (+ `AMDrivingConfig`, `ConstantDrivingConfig`)  
  canonical default external input model used by engine and metadata snapshot.
- `SimulationResult`  
  typed payload returned by the engine and consumed by analyzer/plots.

The goal is to avoid loose dictionary key coupling and make APIs explicit.

## Parameter Layer

`LaNMM_parameters.jl` is the single source of truth for:

- synapse type map
- base connectivity constants
- condition modifiers (`healthy`, `mci`, `ad`, `psychedelics`, `mci+psy`, `ad+psy`)
- default driving model

Current default driving profile:
- `e1`: `multiscale` (baseline swept by `mu_p1_values`)
- `e2`: `constant` (baseline swept by `mu_p2_values`)
- `pv`: `constant` (default baseline `0.0`)

Public functions:

- `get_column_parameters(; condition)`
- `build_synapse_parameter_arrays(column_params)`
- `get_default_driving_model()`

## Engine Layer

`LaNMM_Engine.jl` contains:

- `LaNMMParams` (runtime ODE parameter pack)
- `lanmm_ode!` (14 second-order synapses, 28 states)
- `run_unified_simulation(intrinsic::IntrinsicConfig, job::JobConfig; mu_e1, mu_e2, mu_pv, driving)`
- compatibility wrapper `run_unified_simulation(; ...)`

Numerical robustness:

- primary solver: `Tsit5`
- fallback solver: `Rosenbrock23`
- retcode checks
- discard bounds checks
- output-drive alignment on final simulation time grid

## Analyzer Layer

`LaNMM_Analyzer.jl` is analysis-only:

- `compute_s2e_e2e`
- `analyze_sweep_couplings`
- `compute_band_power`
- `analyze_sweep_power`
- `compute_peak_frequency`
- `analyze_peak_frequencies`
- `compute_peix`
- `sweep_peix`

No plotting side effects are implemented here.

## Plot Layer (Phase 3)

`LaNMM_Plots.jl` contains all plotting:

- `plot_sim_results`
- `plot_coupling_heatmaps`
- `plot_power_heatmaps_bottom_cbar`
- `plot_peix_heatmaps`
- `plot_frequency_heatmaps`

This keeps visualization concerns decoupled from computation.

## Sweep / Orchestration Layer

`LaNMM_Sweep.jl` runs full jobs end-to-end:

- validates config (`validate_configs`)
- writes full run snapshot (`save_parameters_snapshot`)
- runs nominal simulation + diagnostics
- executes parameter sweep with progress bar / ETA / average sec per sim / finish clock
- computes analysis products
- saves plots and serialized outputs

## Run Entry Scripts

Each `run_job_*.jl` file should remain simple:

1. include the Julia pipeline files
2. define `intrinsic_params::IntrinsicConfig`
3. define `job_params::JobConfig`
4. call `run_sweep_job(intrinsic_params, job_params; driving=driving_params)`

## Output Artifacts per Run

The timestamped output directory contains:

- `parameters.txt` (TOML-formatted full parameter snapshot)
- `inputs.png`, `v.png`, `u_external.png`, `psd.png`
- `couplings.png`
- `power.png`
- `peix.png`
- `frequency_heatmaps.png`
- `sweep_results.jls`
- `analysis_results.jls`

## Practical Extension Guide

- Change intrinsic model constants/conditions -> `LaNMM_parameters.jl`
- Change ODE or solver behavior -> `LaNMM_Engine.jl`
- Add a metric -> `LaNMM_Analyzer.jl`
- Add/modify figure style -> `LaNMM_Plots.jl`
- Change batch workflow, progress, saving behavior -> `LaNMM_Sweep.jl`
- Change run setup -> `run_job_P1.jl`, `run_job_P2.jl`, or `run_job_PV.jl`

