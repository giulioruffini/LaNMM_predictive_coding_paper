# Julia code

This folder contains the current Julia implementation used for new runs.

Main entry points:

- `run_job_P1.jl`
- `run_job_P2.jl`
- `run_job_PV.jl`

Run from repository root, for example:

```bash
julia --threads auto julia/run_job_P1.jl
```

Pipeline modules:

- `LaNMM_Types.jl`: typed configs/results.
- `LaNMM_parameters.jl`: intrinsic constants + condition modifiers + default drive.
- `LaNMM_Engine.jl`: ODE + single simulation.
- `LaNMM_Analyzer.jl`: metrics only.
- `LaNMM_Plots.jl`: plotting only.
- `LaNMM_Sweep.jl`: orchestration/progress/saving.
