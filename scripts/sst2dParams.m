%% Shared parameters for 2D Single-Target Tracking scripts
%
% Loaded by estimateSST2d.m / estimateSST2dLfm.m (directly) and by the tune*
% functions via collectParams('sst2dParams'), so they agree on plant, network,
% and estimator settings. Scripts that need script-specific overrides (e.g.
% totalRuns) set them after sourcing this file.

%% Plant
dim = 2;                  % Spatial dimension (drives plant ctor + measurement dim m in sstPlant)
T = 1000;                 % Number of simulation steps
Ts = 0.1;                 % Sampling period
outputNoiseStd = 10;
turnRate = 0;             % 0 = straight line; try 0.005 rad/s for coordinated turn
x0 = [25 25 1000 1000]';  % Initial true state: [vx vy px py]

%% Network
nodeCount = 100;
sensorCount = 100;

% The stochastic trigger implementations can't handle non sensor nodes properly
% sensorCount = 20;

% nodeCount = 20;
% sensorCount = 4;

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
% Z in S^m_{++} per Han et al. 2015 — m=2 for this 2D sensor model)
errorNormWeightsClosed = 1e-4 * eye(2);  % Z (closed-loop, innovation space)
errorNormWeightsOpen = 1e-6 * eye(2);    % Y (open-loop, raw-measurement space)

%% Initial estimate
x0_hat = x0;
P0 = 1e3 * eye(size(x0, 1));  % No prior; swap to 1e-6*I for "perfect knowledge"

%% Monte Carlo
totalRuns = 200;
% totalRuns = 100;
% totalRuns = 10;
% totalRuns = 2;

totalTuneRuns = 50;

%% Tuning grids (consumed by the tune* functions; see tuning/)
% Grids are filter-specific: DKF and RDKF explore different (beta, delta)
% ranges, and the SDKF and SRDKF stochastic triggers live on very different z
% scales, so each gets its own grid rather than one shared sweep.

% DKF event-trigger grid (beta x delta; alpha fixed). See tuneDKF.
tuneDkfBetaGrid   = [0.1, 0.2, 0.5, 0.9];
tuneDkfDeltaGrid  = [0.01, 0.05, 0.1, 0.5, 1.0];
tuneDkfAlphaFixed = 1;

% RDKF event-trigger grid (beta x delta; alpha and robust tol b fixed). See tuneRDKF.
tuneRdkfBetaGrid   = [0.1, 0.5, 1, 2.5];
tuneRdkfDeltaGrid  = [0.1, 0.5, 1, 2.5];
tuneRdkfAlphaFixed = 1;
tuneRdkfBFixed     = 0.05;

% Stochastic-trigger weight scale z (Z = z * eye(m)). SDKF and SRDKF differ by
% orders of magnitude. The SRDKF closed and open triggers get separate grids
% (they weight the innovation vs the raw measurement, so their z scales differ);
% both retain the previously-shared values here so 2D behaviour is unchanged.
% See tuneSDKFClosed / tuneSRDKFClosed / tuneSRDKFOpen.
tuneSdkfZGrid       = [1e-5, 2e-5, 3e-5, 4e-5, 4.5e-5, 5e-5, 6e-5, 7e-5, 1e-4];
tuneSrdkfClosedZGrid = [0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300];
tuneSrdkfOpenZGrid   = [0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300];

% Robust KL-tolerance b sweep on least-favorable data (shared by the *bLfm
% functions). b = 0 is the anchor (robust layer disabled); log grid through and
% just past the shipped 0.05.
tuneBGrid = [0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 5e-2, 1e-1];
