# 3D Single-Target Tracking Extension — Implementation Plan

Status: **implemented (2026-07-19).** Steps 1–6 done and smoke-tested; the full
Monte Carlo (step 7, at the full `totalRuns`) is left to run manually. See the
"Implementation status" note at the end.

Author scope agreed with Mauricio (2026-07-19):

- This is a **3D extension of the existing framework**, not a bare reproduction of
  the paper's figures.
- **SDKF and SRDKF must keep running** (the stochastic-trigger filters stay in the
  comparison set, alongside RDKF / RDKFLOC / SRDKFLOC / DKF / DSEACP / CKF / CRKF).
- **Drop** the coordinated-turn generalization — 3D dynamics are a plain triple
  integrator (constant velocity), no turn-rate coupling.
- Keep **`sensorCount = 100` (= `nodeCount`)** for now. Every node is a sensor, so
  the non-sensor / relay-node code paths stay dormant and the known
  stochastic-trigger relay bug is not exercised.

---

## 1. Background: what already exists

The paper (Ghion & Zorzi, 2023, *Robust distributed Kalman filtering with
event-triggered communication*) evaluates its algorithms in **Section 6 on a 3D
target-tracking model** (Eq. 29). What lives in this repo is a **2D adaptation** of
that model:

- `plants/SingleTarget2dModel.m` — `n = 4` state `[vx vy px py]`, coordinated-turn
  `A`, each sensor measures **one** position coordinate but occupies a **2-row**
  measurement block (one row is identically zero).
- `scripts/estimateSST2d.m`, `scripts/estimateSST2dLfm.m` — Monte Carlo drivers
  (nominal data and least-favorable-model data).
- `scripts/sst2dParams.m` — shared plant / network / estimator parameters.

The heavy machinery is **already dimension-agnostic** and needs **no change**:

- `plants/LeastFavorableModel.m` — derives every size from `plant` (`n`, `p`, `m`).
- `predict/findTheta.m`, `networks/calcFusionWeights.m`, `utils/bisect.m`,
  `utils/calculateRmse.m`.
- The information-form recursions inside every estimator (`updateOmega`,
  `getLocalPriors`, `updateGlobalPriors`, fusion) — all use `self.n` and full-size
  matrices.
- `networks/createSpatialNetwork.m` — places nodes in a 2D plane purely for
  **connectivity**; this is independent of the target's state dimension and stays
  as-is.

## 2. The single structural coupling: the 2-row sensor block

The **only** place the state/measurement dimension is hardcoded is the assumption
that **each sensor occupies a 2-row measurement block**, written as
`idx = (2*i - 1):(2*i)`. Node index `i` is used directly as the sensor's block
index, which is valid because `createSpatialNetwork` makes sensors the first `S`
nodes.

Occurrences to generalize:

| File | Lines | Role |
|---|---|---|
| `estimators/RDKF.m` | 113 | correction-step measurement slice |
| `estimators/DKF.m` | 104 | correction-step measurement slice |
| `estimators/DSEACP.m` | 81 | correction-step measurement slice |
| `estimators/SDKF.m` | 108, 135, 173 | correction, trigger, silent-neighbor recon |
| `estimators/SRDKF.m` | 126, 158, 199 | correction, trigger, silent-neighbor recon |
| `scripts/estimateSST2d.m` | 29 | local-observability check loop |
| `scripts/estimateSST2dLfm.m` | 34 | local-observability check loop |
| `utils/computeLocalTolerances.m` | 88–89 | per-node `b^i` slice (`pI = 2`, row map) |

`RDKFLOC` / `SRDKFLOC` inherit their measurement handling from `RDKF` / `SRDKF`, so
fixing the parents covers them.

Centralized filters (`CKF`, `CRKF`) use the **full** stacked `C` and do no per-node
slicing — no correctness change, only plotting (Section 6).

## 3. Target 3D model (paper Eq. 29, extended to our framework)

State (ordering keeps velocities first, then positions, mirroring the 2D layout):

```
x = [vx vy vz px py pz]'  in R^6          n = 6
```

Continuous triple integrator, discretized at Ts = 0.1. The paper uses the Euler
form directly:

```
Phi = [ 0_3   0_3 ;
        I_3   0_3 ]                       (6x6)
A   = I_6 + Ts * Phi                      (Ts = 0.1  ->  A = I_6 + 0.1 Phi)
B   = sqrt(0.001) * I_6
w_t ~ N(0, I_6)                           => Q = B B' = 0.001 I_6
```

Decision: **match the paper's `A = I + Ts*Phi` and `B = sqrt(0.001) I`** rather than
reusing the 2D plant's exact `c2d` + unit-noise convention. (The 2D plant used
`Bc = eye`, i.e. unit process noise, and exact `c2d`; the paper's numbers differ and
we follow the paper for the 3D model.)

**Measurement model — the important extension.** Each sensor measures **2 of the 3
position coordinates**, laid out as a **3-row block with one zero row** (so
`p_i = 3`):

- horizontal-horizontal: `C^i = [ 0_{3x3}  diag(1,1,0) ]`  (measures px, py)
- horizontal-vertical:   `C^i = [ 0_{3x3}  diag(0,1,1) ]`  (measures py, pz)

Split the sensors between the two types (e.g. first half horizontal-horizontal,
second half horizontal-vertical), analogous to the 2D `C_i_a` / `C_i_b` split.

**Measurement noise (permutation-scrambled, per sensor):**

```
R0 = 0.5 * diag(1, 4, 7)                  (3x3)
R^i = sqrt(k) * P_i R0 P_i'               P_i a random 3x3 permutation per sensor
```

`P_i` is drawn once per sensor at construction and **stored** on the plant for
reproducibility (so a re-simulated plant yields the same `R^i`). Pick `k` to give a
noise scale comparable to the current 2D setup (2D used `outputNoiseStd = 10`,
i.e. `R = 100 I`); `k` is a tunable knob — start with `k = 1` and calibrate against
RMSE sanity ranges. `D^i` = a matrix square root (Cholesky) of `R^i` so that the
plant's additive-noise form `y = C x + D v`, `v ~ N(0,I)` reproduces `R^i`.

`x0 ~ N(0, I_6)` in the paper; we keep our convention of a fixed `x0` plus
`P0 = 1e3 * eye(n)` prior (already generic via `size(x0,1)`), and set a 6-vector
`x0` in params.

## 4. Design decision: uniform per-sensor block width `senBlock`

Introduce a single scalar **`senBlock`** = measurement rows per sensor (`= 2` in 2D,
`= 3` in 3D). This preserves the existing **homogeneous-sensor** assumption (all
sensors have the same `p_i`), which is already implicit in the `2*i-1:2*i` pattern.

- Derive it once per estimator constructor as
  `senBlock = plant.p / sum(G.Nodes.isSensor)` (integer; assert it divides evenly).
- Replace every `idx = (2*i - 1):(2*i)` with a block starting at
  `senBlock*(i-1) + 1` of width `senBlock`:
  `idx = senBlock*(i-1) + (1:senBlock)`. Same substitution for the neighbor slice
  `jdx` in SDKF/SRDKF.
- In `computeLocalTolerances.m`, replace `pI = 2` with `pI = senBlock` and the row
  map `n + (2*i-1), n + 2*i` with `n + senBlock*(i-1) + (1:senBlock)`.

This is deliberately **not** a full heterogeneous-`p_i` refactor (per-node index
maps). Homogeneous `senBlock` is the minimal change that makes the shared estimators
work for both 2D and 3D and keeps the diff small. A heterogeneous refactor is noted
as future work in Section 8.

## 5. New / changed files

### New files

1. **`plants/SingleTarget3dModel.m`**
   - `n = 6`; `A = I_6 + Ts*Phi`; `B = sqrt(0.001) I_6`.
   - Per-sensor 3x6 `C^i` alternating horizontal-horizontal / horizontal-vertical;
     stack into `self.C` (`p = 3*S`).
   - Per-sensor `R^i = sqrt(k) P_i R0 P_i'`; store the drawn `P_i` (and/or the
     assembled block-diagonal `R`, `D`) as properties for reproducibility.
   - `stateEq`, `outputEq`, `simulate` copied from the 2D model (drop `omega`).
   - Expose `p` so estimators can compute `senBlock = p / S`.
   - `plotTrajectory` uses `plot3` over position rows 4:6 (see Section 6).

2. **`scripts/sst3dParams.m`** (clone of `sst2dParams.m` with 3D values)
   - `x0` a 6-vector `[vx vy vz px py pz]'` (e.g. velocities + initial position).
   - `errorNormWeightsClosed = 1e-4 * eye(3)`, `errorNormWeightsOpen = 1e-6 * eye(3)`
     (stochastic-trigger `Z` lives in measurement space, now `p_i = 3`).
   - Drop `turnRate`.
   - `nodeCount = sensorCount = 100`, `maxLength = 5000` (unchanged).
   - `klTolerance = lfmKlTolerance = 0.05`, DKF trigger params unchanged.
   - Add a `k` (measurement-noise scale) parameter for `R^i`.
   - `T`: keep the 2D default for now (1000) or set 250 to match the paper — call
     it out as a knob.

3. **`scripts/estimateSST3d.m`** — near-verbatim clone of `estimateSST2d.m`:
   - source `sst3dParams`;
   - `plant = SingleTarget3dModel(Ts, sensorCount, ...)` (no `turnRate`);
   - everything else (estimator construction, Monte Carlo loop, `saveRun`) is
     dimension-agnostic and unchanged.

4. **`scripts/estimateSST3dLfm.m`** — near-verbatim clone of `estimateSST2dLfm.m`;
   `LeastFavorableModel` needs no change.

### Changed files (the `senBlock` generalization from Section 4)

- `estimators/RDKF.m`, `DKF.m`, `DSEACP.m`, `SDKF.m`, `SRDKF.m`
  (add `senBlock` property + replace hardcoded slices).
- `utils/computeLocalTolerances.m` (`pI` + row map).
- `scripts/estimateSST3d.m` / `estimateSST3dLfm.m` local-observability loop uses
  `senBlock` instead of `2`.

Note: after generalization the 2D scripts still work unchanged (`senBlock` resolves
to 2), so this is a non-breaking refactor.

## 6. Plotting (cosmetic, pervasive)

Trajectory plotting hardcodes 2D position rows `(3,:)` / `(4,:)`. In 3D, position is
rows `4:6`. Affected: `plants/SingleTarget2dModel.m` (n/a, keep), the new 3D plant,
`estimators/{RDKF,DKF,SDKF,SRDKF,DSEACP,CKF,CRKF}.m` `plotTrajectory`,
`scripts/plotSST2dTrajectories.m`, `scripts/plotSST2dRun.m`.

Options:
- Minimal: add 3D-aware trajectory plots (`plot3` over rows 4:6) in the new plant
  and a new `scripts/plotSST3dTrajectories.m` / `plotSST3dRun.m`; leave the 2D
  plotters alone.
- Cleaner: parametrize a `posIdx` (position row indices) so one plotter serves both.

RMSE / transmission-rate plots are dimension-agnostic (they operate on scalar
metrics) and need no change.

## 7. Validation checklist

1. **System checks pass**: `assertStabilizable(A,B)`, `assertDetectable(A,C)`, and
   per-sensor `assertLocallyObservable(A, C_i, i)` for the 3D `C_i` blocks (each
   sensor observes 2 of 3 positions; confirm collective observability holds — pair
   `(A,C)` observable per paper Assumption A3).
2. **`senBlock` divides evenly**: `plant.p == senBlock * sensorCount`.
3. **2D regression**: existing `estimateSST2d` / `estimateSST2dLfm` still run and
   produce unchanged RMSE (the refactor must be behavior-preserving for `senBlock=2`).

> **Keep smoke / regression runs small.** Override `totalRuns = 3` (or similar) for
> every smoke and regression run — the full Monte Carlo count is slow. Set it after
> sourcing the params file (the LFM script already does this via a
> `totalRuns = 5; % TEMP smoke override` line). Only bump to the full count for the
> final validation run in step 7.
4. **LFM sanity**: `LeastFavorableModel(plant3d, P0, 0.05, T)` constructs without
   `findTheta`/`chol` failures; if it errors, reduce `b` or check `Q`/`R` scaling.
5. **`computeLocalTolerances`**: `0 < b^i <= bGlobal` for all sensors (the built-in
   warning guards this); steady-state convergence assertions pass.
6. **Smoke run**: small `totalRuns` (e.g. 3) end-to-end on both 3D scripts, with all
   filters (RDKF, RDKFLOC, DKF, SDKF closed/open, SRDKF closed/open, SRDKFLOC,
   DSEACP, CKF, CRKF) producing finite RMSE.
7. **Full Monte Carlo**: bump `totalRuns` and confirm robust filters (RDKF/SRDKF/LOC)
   beat DKF on LFM data, as in the paper.

## 8. Explicitly out of scope / future work

- Coordinated-turn dynamics in 3D (dropped per scope).
- Relay (non-sensor) nodes with `sensorCount < nodeCount` — the stochastic-trigger
  filters (`SDKF`/`SRDKF`) have a known issue with non-sensor nodes
  (`sst2dParams.m:19`); keeping `sensorCount = 100` avoids it. Revisit if an 80/20
  sensor/relay split (paper's actual topology) is wanted later.
- Heterogeneous per-sensor `p_i` (per-node index maps instead of a uniform
  `senBlock`).
- Calibrating `k` / `R0` scaling to a target RMSE regime.

## 9. Suggested order of work

1. Generalize the `senBlock` slice in estimators + `computeLocalTolerances` (keep 2D
   green — regression test).
2. Write `SingleTarget3dModel.m`.
3. Write `sst3dParams.m`.
4. Clone the two driver scripts.
5. 3D-aware plotting.
6. Smoke run -> validation checklist -> full Monte Carlo.

## 10. Implementation status (2026-07-19)

Done and smoke-tested:

- **`senBlock` generalization** — added a `senBlock` property (derived as
  `plant.p / S`, with a divisibility assert) to `RDKF`, `DKF`, `DSEACP`, `SDKF`,
  `SRDKF`, and replaced every `(2*i-1):(2*i)` / `(2*j-1):(2*j)` slice with
  `senBlock*(i-1)+(1:senBlock)`. `computeLocalTolerances.m` derives `senBlock` the
  same way and uses `pI = senBlock` + the generalized row map. `RDKFLOC` / `SRDKFLOC`
  inherit the change. The substitution is provably identical for `senBlock = 2`
  (`2(i-1)+[1,2] = [2i-1,2i]`), so the 2D path is unchanged.
- **New files**: `plants/SingleTarget3dModel.m` (n=6 triple integrator,
  `A = I + Ts*Phi`, `B = sqrt(0.001) I`, per-sensor 3-row C blocks alternating
  hh/hv, `R^i = sqrt(k) P_i R0 P_i'` with stored permutations),
  `scripts/sst3dParams.m`, `scripts/collectParams3d.m`, `scripts/estimateSST3d.m`,
  `scripts/estimateSST3dLfm.m`, `scripts/plotSST3dRun.m` (delegates to the
  dimension-agnostic `plotSST2dRun`), `scripts/plotSST3dTrajectories.m`
  (`plot3` over rows 4:6).

Smoke results (T=60, sensorCount=100, MC=1): all 12 filters produce finite RMSE
on both nominal and LFM data, `0 < b^i <= bGlobal` holds, and on LFM data the
robust filters beat DKF (RDKF 1.02 vs DKF 1.32) — the paper's expected ordering.
The 2D regression smoke passes with `senBlock = 2`.

Notes / caveats:

- **Local observability**: each 3D sensor measures 2 of 6 states, so
  `assertLocallyObservable` *warns* per sensor (unobservable mode at lambda=1) —
  expected; collective detectability is what the algorithm needs and
  `assertDetectable(A,C)` passes. Same situation as 2D.
- The V/H/L "not at steady state" warnings in the smoke run are an artifact of the
  tiny smoke `T`; they vanish at the full `T`.
- **`noiseScale` (k)** defaults to 1 (low measurement noise); calibrating it to a
  target RMSE regime is still open (Section 8).
- **Still to run manually**: the full Monte Carlo (step 7) at the full `totalRuns`
  via `estimateSST3d` / `estimateSST3dLfm` (set `totalRuns = 3` first for a quick
  end-to-end check).
