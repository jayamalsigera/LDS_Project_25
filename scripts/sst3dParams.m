%% Shared parameters for 3D Single-Target Tracking scripts
%
% Loaded by estimateSST3d.m and estimateSST3dLfm.m (and via collectParams by the
% tune* functions) so they agree on plant, network, and estimator settings.
% Scripts that need script-specific overrides (e.g. totalRuns) set them after
% sourcing this file.

%% Plant
dim = 3;                  % Spatial dimension = measurement rows per sensor

Ts = 0.1;                 % Sampling period

T = 250 / Ts;             % Number of simulation steps (from paper)
% T = 250;                  % Smoke tests

noiseScale = 10;           % Measurement-noise scale k: R^i = sqrt(k) P_i R0 P_i'

% Initial true state [vx vy vz px py pz]. Following the paper, x0 is random with
% zero mean and covariance V_{0|-1} = I, redrawn every Monte Carlo run (the
% estimate* scripts call sampleX0 inside their run loops) so the initial
% condition contributes its share of the MC variance.
stateDim = 2 * dim;
V0 = eye(stateDim);                 % V_{0|-1}
sampleX0 = @() randn(stateDim, 1);

% Single fixed draw, kept only for the tune* functions, which still simulate
% from one deterministic initial state (see tuning/).
x0 = sampleX0();

%% Network

connTarget = 0.04;   % Percentage of node connections (excluding self loops);

nodeCount = 100;

sensorCount = 20;     % Similar scenario to Ghion & Zorzi (2023)
% sensorCount = 100;      % For correct behavior of the Stochastic event estimators

% The stochastic trigger implementations can't handle non sensor nodes properly,
% so keep every node a sensor for now (sensorCount = nodeCount).

maxLength = 5000;         % Spatial extent for random network generation

%% Estimators
consensusSteps = 3;       % DSEACP consensus iterations (L)
% consensusSteps = 1;       % Smoke test

% Event-trigger parameters for DKF/RDKF/SRDKF
dkfAlpha = 10;
dkfBeta = 0.2;
dkfDelta = 0.5;

% KL-divergence tolerance for robust filters (b); LFM data is generated at the
% same radius so the filter defends exactly the mismatch present in the data
% (Ghion & Zorzi 2023, Section 6, uses b = 0.05 for both).
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
% Prior matched to the distribution x0 is drawn from: mean E[x0] = 0 and
% covariance V_{0|-1}. Do not hand the filters the realized x0 -- that is a
% "cheating" prior and it makes the transient meaningless.
x0_hat = zeros(stateDim, 1);
P0 = V0;

%% Monte Carlo
totalRuns = 200;
% totalRuns = 10;   % smoke tests

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
