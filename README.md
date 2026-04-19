# Robust Distributed Kalman Filtering (RDKF) – MATLAB Project

This project implements and compares several distributed state-estimation algorithms for sensor networks, including robust and event-triggered variants of the Distributed Kalman Filter. Monte Carlo experiments target 2D single-target tracking over a spatially distributed sensor network, with optional least-favorable-model (LFM) perturbation for stress-testing robust filters.

Implementation follows:

- *Ghion & Zorzi (2023): Robust Distributed Kalman Filtering with Event-Triggered Communication*
- *Battistelli et al. (2018): A Distributed Kalman Filter with Event-Triggered Communication and Guaranteed Stability*
- *Han et al. (2015): Stochastic Event-Triggered Sensor Schedule for Remote State Estimation*
- *Battistelli & Chisci (2014): Stability of Consensus Extended Kalman Filtering for Distributed State Estimation*
- *Levy & Nikoukhah (2013): Robust State-Space Filtering Under Incremental Model Perturbations* (least-favorable model, Section V)

---

## Project Structure

- `estimators/` — Filter implementations: CKF (centralized baseline), DKF, DSEACP, RDKF, and SRDKF.
- `plants/` — Target-dynamics models: the nominal 2D single-target plant and the least-favorable-model data generator, plus stabilizability/detectability checks.
- `networks/` — Spatial sensor-network generation, consensus-weight computation, connectivity diagnostics, and graph visualization.
- `predict/` — Robust prediction primitives used by RDKF/SRDKF (risk-sensitivity bisection, robust fusion, no-transmit prediction).
- `fusion/` — Stochastic- and deterministic-trigger fusion conditions and data fusion.
- `utils/` — General helpers (RMSE, bisection, Loewner-order checks, `parfor` progress bar).
- `scripts/` — Top-level experiments (`estimateSST2d`, `estimateSST2dLfm`, `tuneDKF`, `tuneRDKF`), the shared parameter file `sst2dParams.m`, and the save/load/plot tooling that makes Monte Carlo runs reproducible.
- `papers/` — Reference material.
- `legacy/`, `idkf/`, `ssm/` — Archived or in-progress code outside the main pipeline.

---

## Running experiments

### 1. Open the MATLAB project

Double-click `lds-proj.prj` (this puts all subfolders on the path).

### 2. Configure shared parameters

All four experiment scripts source `scripts/sst2dParams.m`. Edit that file to change anything shared across experiments.

Script-specific settings (e.g. `totalRuns`, hyperparameter grids) stay inside each script, immediately after the `sst2dParams;` call — override shared values there when needed.

### 3. Run a script

From the MATLAB Command Window (any of):

```matlab
estimateSST2d        % nominal-plant Monte Carlo
estimateSST2dLfm     % least-favorable-model Monte Carlo
tuneDKF              % DKF event-trigger sweep
tuneRDKF             % RDKF event-trigger + KL-tolerance sweep
```

Each script:

1. Builds a spatial sensor network and a target plant (and, for the LFM variant, precomputes the worst-case data generator).
2. Runs the configured Monte Carlo trials in parallel via `parfor`.
3. **Saves the run** to `results/<scriptName>_T<T>_N<nodes>s<sensors>_b<klTol>_runs<N>_<timestamp>.mat`.
4. Reloads the saved run and plots from it (so fresh and reloaded runs produce identical figures).

---

## Saving and reloading results

Every experiment run is persisted automatically. The filename encodes the most commonly varied parameters; the full parameter set is stored inside the `.mat`.

Example: `results/estimateSST2dLfm_T1000_N100s20_b0.05_runs100_20260419-153045.mat`

The file contains a single `runData` struct with:

| Field        | Contents                                                               |
|--------------|------------------------------------------------------------------------|
| `script`     | Name of the producing script                                            |
| `timestamp`  | ISO-ish local time (`yyyyMMdd-HHmmss`)                                  |
| `gitSha`     | Short git SHA at save time (best-effort)                                |
| `params`     | Every variable defined in `sst2dParams.m`                               |
| `extras`     | Script-specific config (`totalRuns`, hyperparameter grids/baselines)    |
| `netGraph`   | Sensor/communication graph                                              |
| `results`    | Per-run RMSE and TX-rate arrays (and sweep tables for tune scripts)     |
| `samples`    | Showcase sample objects used by trajectory plots (estimate scripts)     |

### Reload and replay figures

```matlab
plotSST2dRun('results/estimateSST2dLfm_...mat');   % estimate-script runs
% or
plotTuneRun('results/tuneRDKF_...mat');            % tune-script runs
```

`plotSST2dRun` automatically detects LFM runs (via the presence of `samples.nomSample`) and draws the nominal-vs-LFM trajectory overlay. `plotTuneRun` detects RDKF runs (via a `b` field in `configs`) and emits the additional KL-tolerance panel.

---

## Output

Figures produced by the estimate scripts:

- Network layout
- True trajectory (LFM scripts additionally overlay the nominal trajectory for contrast)
- Per-estimator estimated trajectory vs. truth
- RMSE vs. time (log scale), averaged over Monte Carlo runs
- Transmission rate vs. time, averaged over Monte Carlo runs

Figures produced by the tune scripts:

- RMSE and TX-rate vs. time for each swept parameter (one subplot pair per parameter)
- RMSE-vs-TX-rate tradeoff scatter, color-coded by which parameter is being swept

Each estimator sample object also carries:

```matlab
sample.RMSE     % RMSE over time
sample.txRate   % Transmission rate over time (DKF/RDKF/SRDKF)
sample.X_hat    % Estimated state at each node (distributed filters)
sample.x_hat    % Estimated state (centralized filter, CKF)
```
