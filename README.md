# Robust Distributed Kalman Filtering (RDKF) – MATLAB Project

This project implements and compares several distributed state estimation algorithms for sensor networks, including robust and event‑triggered variants of the Distributed Kalman Filter (DKF). The implementation is based on:

- *Ghion & Zorzi (2023): Robust Distributed Kalman Filtering with Event‑Triggered Communication*
- *Battistelli, et al., (2018): A distributed Kalman filter with event-triggered communication and guaranteed stability*
- *Han, et al. (2015): Stochastic Event-Triggered Sensor Schedule for Remote State Estimation*
- *Batistelli & Chisci (2014): Stability of Consensus Extended Kalman Filtering for Distributed State Estimation*

The implementation is done for the problem of 2D single‑target tracking with noisy measurements and a spatially distributed sensor network.

---

## Project Structure

### `estimators/`
MATLAB classes implementing different filtering algorithms:

- **CKF.m** – Centralized Kalman Filter
- **DKF.m** – Distributed Kalman Filter (Battistelli, et al.)
- **DSEACP.m** – Distributed State Estimation with Average Consensus Protocol (Batistelli & Chisci)
- **RDKF.m** – Robust DKF with event‑triggered communication (Ghion & Zorzi)
- **SRDKF.m** – Stochastic Robust DKF (stochastic triggering variant)

```

---

### `networks/`
Utilities for generating and analyzing communication graphs:

- **createSpatialNetwork.m** – Random geometric graph with designated sensor nodes
- **calcMetropolisWeights.m** – Metropolis–Hastings consensus weights
- **calcFusionWeights.m** – Uniform consensus weights (as used in the paper)
- **plotNetwork.m** – Visualization of node positions and edges
- **getConnectivityPercentage.m** – Connectivity diagnostics

---

### `plants/`
Plant models used for simulation:

- **SingleTarget2dModel.m** – 2D kinematic model with position/velocity state and noisy sensor outputs

This class provides:
- system matrices (A, B, C, D)
- simulation of true trajectory and measurements

---

### `scripts/`
Top‑level scripts for running experiments:

- **estimateSST2d.m** – Main simulation script (Monte Carlo runs, RMSE plots, Tx‑rate plots)
- **testConnectivityPercentage.m** – Utility for checking network connectivity statistics

---

### `utils/`
General helper functions:

- **loewnerBetweenEig.m** – Checks Loewner ordering between matrices (used in event triggers)

---

### Other files
- **lds-proj.prj** – MATLAB project file
- **.vscode/** – Editor configuration

---

## To run project:

### 1. Open the MATLAB project

- Double‑click `lds-proj.prj

---

### 2. Run the main simulation script

From MATLAB:

```matlab
scripts/estimateSST2d
```

This script:

1. Creates a spatial sensor network
2. Simulates the 2D target trajectory
3. Runs CKF, DKF, DSEACP, RDKF, and SRDKF
4. Computes RMSE and transmission rates
5. Plots:
   - true vs. estimated trajectories
   - RMSE over time
   - Tx‑rate over time

---

### 3. Adjustment of parameters

Inside `estimateSST2d.m`, the following parameters may be adjusted

- **Network size**
  ```matlab
  nodeCount = 20;
  sensorCount = 4;
  ```

- **Event‑trigger parameters**
  ```matlab
  dkfAlpha = 10;
  dkfBeta  = 0.2;
  dkfDelta = 0.5;
  ```

- **Robustness tolerance**
  ```matlab
  rdkfb = 0.05;
  ```

- **Initial covariance**
  ```matlab
  P0 = 1e3 * eye(4);
  ```

- **Number of Monte Carlo runs**
  ```matlab
  totalRuns = 10;
  ```

---

## Output

The following outputs are plottted by estimateSST2d

- **Network plot**
- **True and estimated trajectories**
- **RMSE comparison (log scale)**
- **Transmission rate comparison**

Each estimator object also stores:

```matlab
result.RMSE     % RMSE over time
result.txRate   % Transmission rate over time (for DKF/RDKF/SRDKF)
result.X_hat    % Estimated states
```

##Issues and todos

As this is only a preliminary version of the project, several functions are not implemented:
- Several checks are not present in the script: Connectivity, Stabilizability and Detetability
- The plant used to simulate the dynamics for all the plants is a nominal model, ie. not stochastically generated or a worst case model. The same goes for initial conditions.
- The matrix Psi in the robust distributed Kalman filter and robust distributed Kalman filter with stochastic trigger becomes ill-conditioned. It has not been possible to determine a mistake in the code compared to the algorithm, so it is possible that the error is due to the initial conditions or how the robust filter works with a nominal plant.
- The way the robust kalman filter is constructed right now, its trigger is always active. Again the source of this has not been found, but the ill-conditioned Psi-matrix would probably mean that a high innovation is seen each iteration leading to exceeding the tolerance.
