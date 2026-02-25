# Python Implementation README

This document mirrors the Julia implementation guide, but for the legacy Python
LaNMM pipeline in `python/`.

## Scope

The Python code is kept for reproducibility and comparison with the Julia
implementation. The active files are:

- `python/lanmmv11.py` (model + simulation + sweep engine)
- `python/lanmm_analyzer.py` (analysis + plotting utilities)
- `python/lanmm_helpers.py` (job orchestration helpers)

## Setup

Create and activate a virtual environment, then install dependencies:

```bash
cd /path/to/LaNMM_predictive_coding_paper
python3 -m venv .venv
source .venv/bin/activate
pip install -r python/requirements.txt
```

Windows (PowerShell):

```powershell
cd C:\path\to\LaNMM_predictive_coding_paper
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -r python/requirements.txt
```

## How to run

Run from inside the `python/` folder so local imports resolve cleanly:

```bash
cd python
python lanmmv11.py
```

`lanmm_helpers.py` can be used to orchestrate sweep jobs and export plots/data
similarly to the Julia workflow.

## Architecture map (Python)

### `lanmmv11.py`

- intrinsic/disease-condition parameter builders
- drive/noise generators (constant, multiscale, AM, pulse)
- core ODE right-hand side (`lanmm_ode_unified`)
- single simulation wrapper (`run_unified_simulation`)
- parameter sweep (`run_parameter_sweep`)

### `lanmm_analyzer.py`

- coupling metrics (S2E/E2E)
- band-power analysis
- PEIX analysis
- peak-frequency maps
- all plotting utilities

### `lanmm_helpers.py`

- higher-level `run_sweep_job(...)` helper
- parameter snapshot serialization
- end-to-end output generation

## Outputs

Typical sweep outputs include:

- `parameters.json`
- `inputs.png`, `v.png`, `u_external.png`, `psd.png`
- `couplings.png`, `power.png`, `peix.png`, `frequency_heatmaps.png`
- `sweep_results.pkl`, `analysis_results.pkl`

## Notes

- Python implementation is monolithic and less modular than Julia.
- If you launch from repo root, use `python -m` style or update `PYTHONPATH`.
- For large sweeps, reduce sweep grid and/or simulation length first.
