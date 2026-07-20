%% Shared parameters for 3D Single-Target Tracking scripts
%
% Loaded by estimateSST3d.m and estimateSST3dLfm.m so they agree on plant,
% network, and estimator settings. Mirrors sst2dParams.m; see the 3D extension
% plan (sst3d-extension-plan.md) for the rationale behind the 3D-specific values.
% Scripts that need script-specific overrides (e.g. totalRuns) set them after
% sourcing this file.

%% Plant
dim = 3;                  % Spatial dimension (drives plant ctor + measurement dim m in sstPlant)
T = 250;                  % Number of simulation steps (paper uses 250; knob)
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
% Z in S^m_{++} per Han et al. 2015 — m=3 for this 3D sensor model)
errorNormWeightsClosed = 1e-4 * eye(3);  % Z (closed-loop, innovation space)
errorNormWeightsOpen = 1e-6 * eye(3);    % Y (open-loop, raw-measurement space)

%% Initial estimate
x0_hat = x0;
P0 = 1e3 * eye(size(x0, 1));  % No prior; swap to 1e-6*I for "perfect knowledge"

%% Monte Carlo
totalRuns = 200;
% totalRuns = 3;   % smoke / regression override (see sst3d-extension-plan.md)

totalTuneRuns = 50;

%% Tuning grids (consumed by the tune* functions; see tuning/)
% Mirrors the sst2dParams tuning section. Grids are filter-specific (DKF vs
% RDKF (beta, delta) ranges differ; SDKF and SRDKF z scales differ by orders of
% magnitude). Starting points for 3D -- the (beta, delta) grids bracket the 3D
% baseline (dkfAlpha = 10, dkfDelta = 0.5); refine after the first sweeps.

% DKF event-trigger grid (beta x delta; alpha fixed). See tuneDKF.
tuneDkfBetaGrid   = [0.1, 0.2, 0.5, 0.9];
tuneDkfDeltaGrid  = [0.05, 0.1, 0.5, 1.0, 2.0];
tuneDkfAlphaFixed = 10;

% RDKF event-trigger grid (beta x delta; alpha and robust tol b fixed). See tuneRDKF.
tuneRdkfBetaGrid   = [0.1, 0.5, 1, 2.5];
tuneRdkfDeltaGrid  = [0.1, 0.5, 1, 2.5];
tuneRdkfAlphaFixed = 10;
tuneRdkfBFixed     = 0.05;

% Stochastic-trigger weight scale z (Z = z * eye(m)). SDKF and SRDKF differ by
% orders of magnitude. See tuneSDKFClosed / tuneSRDKFClosed / tuneSRDKFOpen.
tuneSdkfZGrid  = [1e-5, 2e-5, 3e-5, 4e-5, 4.5e-5, 5e-5, 6e-5, 7e-5, 1e-4];
tuneSrdkfZGrid = [0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300];

% Robust KL-tolerance b sweep on least-favorable data (shared by the *bLfm
% functions). b = 0 is the anchor (robust layer disabled); log grid through and
% just past the shipped 0.05.
tuneBGrid = [0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 5e-2, 1e-1];
