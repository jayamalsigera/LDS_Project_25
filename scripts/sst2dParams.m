%% Shared parameters for 2D Single-Target Tracking scripts
%
% Loaded by estimateSST2d.m, estimateSST2dLfm.m, tuneDKF.m, and tuneRDKF.m
% so they agree on plant, network, and estimator settings. Scripts that need
% script-specific overrides (e.g. totalRuns, hyperparameter grids) set them
% after sourcing this file.

%% Plant
T = 1000;                 % Number of simulation steps
Ts = 0.1;                 % Sampling period
outputNoiseStd = 10;
turnRate = 0;             % 0 = straight line; try 0.005 rad/s for coordinated turn
x0 = [25 25 1000 1000]';  % Initial true state: [vx vy px py]

%% Network
nodeCount = 100;
sensorCount = 20;
% nodeCount = 20;
% sensorCount = 4;

maxLength = 5000;         % Spatial extent for random network generation

%% Estimators
consensusSteps = 3;       % DSEACP consensus iterations (L)

% Event-trigger parameters for DKF/RDKF/SRDKF
dkfAlpha = 1;
dkfBeta = 1;
dkfDelta = 1;

% KL-divergence tolerance for robust filters (b); LFM uses the same b.
klTolerance = 0.05;
lfmKlTolerance = 0.01;

% Stochastic trigger error-norm weight matrices
errorNormWeightsClosed = 0.01 * eye(length(x0));  % Z from the paper
errorNormWeightsOpen = 0.01 * eye(2);             % Y from the paper

%% Initial estimate
x0_hat = x0;
P0 = 1e3 * eye(size(x0, 1));  % No prior; swap to 1e-6*I for "perfect knowledge"

%% Monte Carlo
totalRuns = 200;
% totalRuns = 100;
% totalRuns = 10;
% totalRuns = 2;
