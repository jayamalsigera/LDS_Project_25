%% Shared parameters for 3D Single-Target Tracking scripts
%
% Loaded by estimateSST3d.m and estimateSST3dLfm.m (and via collectParams by the
% tune* functions) so they agree on plant, network, and estimator settings.
% Scripts that need script-specific overrides (e.g. totalRuns) set them after
% sourcing this file.

%% Plant
dim = 3;                  % Spatial dimension = measurement rows per sensor
T = 1000;                 % Number of simulation steps (paper uses 250; knob)
Ts = 0.1;                 % Sampling period
noiseScale = 1;           % Measurement-noise scale k: R^i = sqrt(k) P_i R0 P_i'
x0 = [25 25 25 1000 1000 1000]';  % Initial true state: [vx vy vz px py pz]

%% Network
nodeCount = 100;
sensorCount = 100;

% The stochastic trigger implementations can't handle non sensor nodes properly,
% so keep every node a sensor for now (sensorCount = nodeCount).

maxLength = 5000;         % Spatial extent for random network generation

%% Estimators
consensusSteps = 3;       % DSEACP consensus iterations (L)

% Event-trigger parameters for DKF/RDKF/SRDKF
dkfAlpha = 10;
dkfBeta = 0.2;
dkfDelta = 0.5;

% KL-divergence tolerance for robust filters (b); LFM data is generated at the
% same radius so the filter defends exactly the mismatch present in the data
% (matches Ghion & Zorzi 2023, Section 6, which uses b = 0.05 for both).
klTolerance = 0.05;
lfmKlTolerance = 0.05;

% Stochastic trigger error-norm weight matrices (both in measurement space,
% Z in S^m_{++} per Han et al. 2015 — m=3 for this 3D sensor model). The
% closed-loop weight is set from the tuneSDKFClosed/tuneSRDKFClosed sweeps:
% z=0.15 puts both closed triggers at ~34% TX, near the RMSE knee (the old
% 1e-4 gave ~0.5% TX, starving the filter — see the recalibrated z grids
% below). The open-loop weight already sits near saturation (~94% TX) at 1e-6.
errorNormWeightsClosed = 0.15 * eye(3);  % Z (closed-loop, innovation space)
errorNormWeightsOpen = 1e-6 * eye(3);    % Y (open-loop, raw-measurement space)

%% Initial estimate
x0_hat = x0;
P0 = 1e3 * eye(size(x0, 1));  % No prior; swap to 1e-6*I for "perfect knowledge"

%% Monte Carlo
totalRuns = 200;
% totalRuns = 3;   % smoke / regression override (see sst3d-extension-plan.md)

totalTuneRuns = 50;

%% Tuning grids (consumed by the tune* functions; see tuning/)
% Grids are filter-specific (DKF vs RDKF (beta, delta) ranges differ; SDKF and
% SRDKF z scales differ by orders of magnitude). The (beta, delta) grids bracket
% the baseline (dkfAlpha = 10, dkfDelta = 0.5); refine after the first sweeps.

% DKF event-trigger grid (beta x delta; alpha fixed). See tuneDKF.
tuneDkfBetaGrid   = [0.1, 0.2, 0.5, 0.9];
tuneDkfDeltaGrid  = [0.05, 0.1, 0.5, 1.0, 2.0];
tuneDkfAlphaFixed = 10;

% RDKF event-trigger grid (beta x delta; alpha and robust tol b fixed). See tuneRDKF.
tuneRdkfBetaGrid   = [0.1, 0.5, 1, 2.5];
tuneRdkfDeltaGrid  = [0.1, 0.5, 1, 2.5];
tuneRdkfAlphaFixed = 10;
tuneRdkfBFixed     = 0.05;

% Stochastic-trigger weight scale z (Z = z * eye(m)). The three z-sweeps need
% grids that differ by orders of magnitude because the trigger fires with
% probability 1 - exp(-1/2 z'Zz) evaluated on very different quantities:
% closed-loop weights the innovation y - C x_bar (O(10)), open-loop weights the
% raw measurement y (O(1e3), dominated by position). Grids picked to span
% TX ~2%-95% with the RMSE knee resolved; SRDKF-Closed starts at 0.02 because
% the robust layer destabilises when transmissions are starved below that.
tuneSdkfZGrid       = [3e-3, 8e-3, 2e-2, 4e-2, 8e-2, 0.15, 0.3, 0.6, 1.2];  % SDKF closed
tuneSrdkfClosedZGrid = [2e-2, 4e-2, 8e-2, 0.15, 0.3, 0.6, 1.2, 2.5, 5];      % SRDKF closed
tuneSrdkfOpenZGrid   = [1e-8, 3e-8, 6e-8, 1e-7, 2e-7, 4e-7, 7e-7, 1.5e-6, 3e-6]; % SRDKF open (raw-measurement space)

% Robust KL-tolerance b sweep on least-favorable data (shared by the *bLfm
% functions). b = 0 is the anchor (robust layer disabled); log grid through and
% just past the shipped 0.05.
tuneBGrid = [0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 5e-2, 1e-1];
