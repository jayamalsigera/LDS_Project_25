
function savedPath = estimateSST3d(varargin)
%% Simulations of the 3D Single-Target Tracking Plant

close all;
clc;
rng(42);


%% Parameters
sst3dParams;
P = collectParams(varargin{:});
fields = fieldnames(P);
for f = 1:numel(fields)
  eval([fields{f} ' = P.' fields{f} ';']);
end
totalRuns = P.totalRuns;
sampleX0 = P.sampleX0;
x0_hat = P.x0_hat;
P0 = P.P0;

%% Network Definition
disp("Creating Network")

% netGraph = createSpatialNetwork(nodeCount, sensorCount, maxLength);
netGraph = createRandomNetwork(nodeCount, sensorCount, connTarget, maxLength);

assertConnected(netGraph);

%% Model Simulation
disp("Simulating target dynamics")
plant = SingleTarget3dModel(Ts, sensorCount, noiseScale, T);

assertStabilizable(plant.A, plant.B);
assertDetectable(plant.A, plant.C);

%% System checks (stabilizability, detectability, local observability)
disp("Checking local observability of sensor nodes")
assertLocallyObservable(plant.A, plant.C, sensorCount);

%% Estimators
disp("Initializing estimators")

% Initiate all filters
ckf        = CKF(plant, Ts, T);
crkf       = CRKF(plant, Ts, T, klTolerance);
dseacp     = DSEACP(plant, Ts, T, netGraph, consensusSteps);
dkf        = DKF(plant, Ts, T, netGraph, dkfAlpha, dkfBeta, dkfDelta);
rdkf       = RDKF(plant, Ts, T, netGraph, dkfAlpha, dkfBeta, dkfDelta, klTolerance);
sdkfClosed = SDKF(plant, Ts, T, netGraph, dkfDelta, errorNormWeightsClosed, 'closed');
sdkfOpen   = SDKF(plant, Ts, T, netGraph, dkfDelta, errorNormWeightsOpen,   'open');
srdkfClosed = SRDKF(plant, Ts, T, netGraph, dkfAlpha, dkfBeta, dkfDelta, ...
                    klTolerance, errorNormWeightsClosed, 'closed');
srdkfOpen = SRDKF(plant, Ts, T, netGraph, dkfAlpha, dkfBeta, dkfDelta, ...
                  klTolerance, errorNormWeightsOpen, 'open');

% Local-tolerance robust filters (Algorithm 2). Each constructor computes its
% per-node b^i from the global model of radius klTolerance; this runs serially,
% once. On nominal data these should track their b=0 counterparts (robustness
% earns nothing without model mismatch); the payoff shows on LFM data.
rdkfloc        = RDKFLOC(plant, Ts, T, netGraph, dkfAlpha, dkfBeta, dkfDelta, ...
                         klTolerance, P0);
srdkflocClosed = SRDKFLOC(plant, Ts, T, netGraph, dkfAlpha, dkfBeta, dkfDelta, ...
                          klTolerance, errorNormWeightsClosed, 'closed', P0);
srdkflocOpen   = SRDKFLOC(plant, Ts, T, netGraph, dkfAlpha, dkfBeta, dkfDelta, ...
                          klTolerance, errorNormWeightsOpen, 'open', P0);

bLocal = rdkfloc.b;                        % per-node b^i (same for all LOC filters)
bSens  = bLocal(netGraph.Nodes.isSensor);
fprintf(['Local tolerances b^i (%d sensor nodes): ' ...
         'min=%.3g  median=%.3g  max=%.3g   (global b=%.3g)\n'], ...
        numel(bSens), min(bSens), median(bSens), max(bSens), klTolerance);

%% Monte Carlo simulations
disp("Running Monte Carlo simulations")

% Preallocate RMSE and Transmission rate logs
ckfRmse        = zeros(totalRuns, T + 1);
crkfRmse       = zeros(totalRuns, T + 1);
dseacpRmse     = zeros(totalRuns, T + 1);
dkfRmse        = zeros(totalRuns, T + 1);
rdkfRmse       = zeros(totalRuns, T + 1);
sdkfClosedRmse = zeros(totalRuns, T + 1);
sdkfOpenRmse   = zeros(totalRuns, T + 1);
srdkfClosedRmse = zeros(totalRuns, T + 1);
srdkfOpenRmse   = zeros(totalRuns, T + 1);
rdkflocRmse        = zeros(totalRuns, T + 1);
srdkflocClosedRmse = zeros(totalRuns, T + 1);
srdkflocOpenRmse   = zeros(totalRuns, T + 1);

dkfTxRate        = zeros(totalRuns, T + 1);
rdkfTxRate       = zeros(totalRuns, T + 1);
sdkfClosedTxRate = zeros(totalRuns, T + 1);
sdkfOpenTxRate   = zeros(totalRuns, T + 1);
srdkfClosedTxRate = zeros(totalRuns, T + 1);
srdkfOpenTxRate   = zeros(totalRuns, T + 1);
rdkflocTxRate        = zeros(totalRuns, T + 1);
srdkflocClosedTxRate = zeros(totalRuns, T + 1);
srdkflocOpenTxRate   = zeros(totalRuns, T + 1);

crkfTheta            = zeros(totalRuns, T + 1);
rdkfTheta            = zeros(totalRuns, nodeCount, T + 1);
rdkfThetaBar         = zeros(totalRuns, nodeCount, T + 1);
rdkflocTheta         = zeros(totalRuns, nodeCount, T + 1);
rdkflocThetaBar      = zeros(totalRuns, nodeCount, T + 1);
srdkfClosedTheta     = zeros(totalRuns, nodeCount, T + 1);
srdkfClosedThetaBar  = zeros(totalRuns, nodeCount, T + 1);
srdkfOpenTheta       = zeros(totalRuns, nodeCount, T + 1);
srdkfOpenThetaBar    = zeros(totalRuns, nodeCount, T + 1);
srdkflocClosedTheta  = zeros(totalRuns, nodeCount, T + 1);
srdkflocClosedThetaBar = zeros(totalRuns, nodeCount, T + 1);
srdkflocOpenTheta    = zeros(totalRuns, nodeCount, T + 1);
srdkflocOpenThetaBar = zeros(totalRuns, nodeCount, T + 1);

tic

parfor run = 1:totalRuns
  fprintf('Running simulation %d/%d\n', run, totalRuns);
  s = plant.simulate(sampleX0());

  ckfRmse(run, :)    = ckf.estimate(x0_hat, P0, s.X, s.Y).RMSE;
  crkfRun = crkf.estimate(x0_hat, P0, s.X, s.Y);
  crkfRmse(run, :) = crkfRun.RMSE;
  crkfTheta(run, :) = crkfRun.theta_hist;
  dseacpRmse(run, :) = dseacp.estimate(x0_hat, P0, s.X, s.Y).RMSE;

  dkfRun = dkf.estimate(x0_hat, P0, s.X, s.Y);
  dkfRmse(run, :)   = dkfRun.RMSE;
  dkfTxRate(run, :) = dkfRun.txRate;

  rdkfRun = rdkf.estimate(x0_hat, P0, s.X, s.Y);
  rdkfRmse(run, :)   = rdkfRun.RMSE;
  rdkfTxRate(run, :) = rdkfRun.txRate;
  rdkfTheta(run, :, :) = rdkfRun.theta_hist;
  rdkfThetaBar(run, :, :) = rdkfRun.theta_bar_hist;

  sdkfClosedRun = sdkfClosed.estimate(x0_hat, P0, s.X, s.Y);
  sdkfClosedRmse(run, :)   = sdkfClosedRun.RMSE;
  sdkfClosedTxRate(run, :) = sdkfClosedRun.txRate;

  sdkfOpenRun = sdkfOpen.estimate(x0_hat, P0, s.X, s.Y);
  sdkfOpenRmse(run, :)   = sdkfOpenRun.RMSE;
  sdkfOpenTxRate(run, :) = sdkfOpenRun.txRate;

  srdkfClosedRun = srdkfClosed.estimate(x0_hat, P0, s.X, s.Y);
  srdkfClosedRmse(run, :)   = srdkfClosedRun.RMSE;
  srdkfClosedTxRate(run, :) = srdkfClosedRun.txRate;
  srdkfClosedTheta(run, :, :) = srdkfClosedRun.theta_hist;
  srdkfClosedThetaBar(run, :, :) = srdkfClosedRun.theta_bar_hist;

  srdkfOpenRun = srdkfOpen.estimate(x0_hat, P0, s.X, s.Y);
  srdkfOpenRmse(run, :)   = srdkfOpenRun.RMSE;
  srdkfOpenTxRate(run, :) = srdkfOpenRun.txRate;
  srdkfOpenTheta(run, :, :) = srdkfOpenRun.theta_hist;
  srdkfOpenThetaBar(run, :, :) = srdkfOpenRun.theta_bar_hist;

  rdkflocRun = rdkfloc.estimate(x0_hat, P0, s.X, s.Y);
  rdkflocRmse(run, :)   = rdkflocRun.RMSE;
  rdkflocTxRate(run, :) = rdkflocRun.txRate;
  rdkflocTheta(run, :, :) = rdkflocRun.theta_hist;
  rdkflocThetaBar(run, :, :) = rdkflocRun.theta_bar_hist;

  srdkflocClosedRun = srdkflocClosed.estimate(x0_hat, P0, s.X, s.Y);
  srdkflocClosedRmse(run, :)   = srdkflocClosedRun.RMSE;
  srdkflocClosedTxRate(run, :) = srdkflocClosedRun.txRate;
  srdkflocClosedTheta(run, :, :) = srdkflocClosedRun.theta_hist;
  srdkflocClosedThetaBar(run, :, :) = srdkflocClosedRun.theta_bar_hist;

  srdkflocOpenRun = srdkflocOpen.estimate(x0_hat, P0, s.X, s.Y);
  srdkflocOpenRmse(run, :)   = srdkflocOpenRun.RMSE;
  srdkflocOpenTxRate(run, :) = srdkflocOpenRun.txRate;
  srdkflocOpenTheta(run, :, :) = srdkflocOpenRun.theta_hist;
  srdkflocOpenThetaBar(run, :, :) = srdkflocOpenRun.theta_bar_hist;
end
fprintf('Elapsed: %.3f s\n', toc);

%% Save run
results = struct( ...
  'ckfRmse', ckfRmse, 'crkfRmse', crkfRmse, 'dseacpRmse', dseacpRmse, ...
  'dkfRmse', dkfRmse, 'rdkfRmse', rdkfRmse, ...
  'sdkfClosedRmse', sdkfClosedRmse, 'sdkfOpenRmse', sdkfOpenRmse, ...
  'srdkfClosedRmse', srdkfClosedRmse, 'srdkfOpenRmse', srdkfOpenRmse, ...
  'rdkflocRmse', rdkflocRmse, ...
  'srdkflocClosedRmse', srdkflocClosedRmse, 'srdkflocOpenRmse', srdkflocOpenRmse, ...
  'dkfTxRate', dkfTxRate, 'rdkfTxRate', rdkfTxRate, ...
  'sdkfClosedTxRate', sdkfClosedTxRate, 'sdkfOpenTxRate', sdkfOpenTxRate, ...
  'srdkfClosedTxRate', srdkfClosedTxRate, 'srdkfOpenTxRate', srdkfOpenTxRate, ...
  'rdkflocTxRate', rdkflocTxRate, ...
  'srdkflocClosedTxRate', srdkflocClosedTxRate, 'srdkflocOpenTxRate', srdkflocOpenTxRate, ...
  'crkfTheta', crkfTheta, 'rdkfTheta', rdkfTheta, 'rdkflocTheta', rdkflocTheta, ...
  'srdkfClosedTheta', srdkfClosedTheta, 'srdkfOpenTheta', srdkfOpenTheta, ...
  'srdkflocClosedTheta', srdkflocClosedTheta, 'srdkflocOpenTheta', srdkflocOpenTheta, ...
  'rdkfThetaBar', rdkfThetaBar, 'rdkflocThetaBar', rdkflocThetaBar, ...
  'srdkfClosedThetaBar', srdkfClosedThetaBar, 'srdkfOpenThetaBar', srdkfOpenThetaBar, ...
  'srdkflocClosedThetaBar', srdkflocClosedThetaBar, 'srdkflocOpenThetaBar', srdkflocOpenThetaBar);
extras = struct('totalRuns', totalRuns, 'bLocal', bLocal);
savedPath = saveRun(mfilename, P, extras, netGraph, results);
