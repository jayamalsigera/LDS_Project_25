%% Shared parameters for 3D Single-Target Tracking scripts
%
% Loaded by estimateSST3d.m and estimateSST3dLfm.m so they agree on plant,
% network, and estimator settings. Mirrors sst2dParams.m; see the 3D extension
% plan (sst3d-extension-plan.md) for the rationale behind the 3D-specific values.
% Scripts that need script-specific overrides (e.g. totalRuns) set them after
% sourcing this file.

%% Plant
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
dkfAlpha = 1;
dkfBeta = 0.2;
dkfDelta = 0.1;

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
