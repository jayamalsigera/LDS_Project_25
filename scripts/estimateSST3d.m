
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
sdkf       = SDKF(plant, Ts, T, netGraph, errorNormWeights);
srdkf      = SRDKF(plant, Ts, T, netGraph, klTolerance, errorNormWeights);

% Local-tolerance robust filters (Algorithm 2). Each constructor computes its
% per-node b^i from the global model of radius klTolerance; this runs serially,
% once. On nominal data these should track their b=0 counterparts (robustness
% earns nothing without model mismatch); the payoff shows on LFM data.
rdkfloc    = RDKFLOC(plant, Ts, T, netGraph, dkfAlpha, dkfBeta, dkfDelta, ...
                     klTolerance, P0);
srdkfloc   = SRDKFLOC(plant, Ts, T, netGraph, klTolerance, errorNormWeights, P0);

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
sdkfRmse       = zeros(totalRuns, T + 1);
srdkfRmse      = zeros(totalRuns, T + 1);
rdkflocRmse    = zeros(totalRuns, T + 1);
srdkflocRmse   = zeros(totalRuns, T + 1);

dkfTxRate        = zeros(totalRuns, T + 1);
rdkfTxRate       = zeros(totalRuns, T + 1);
sdkfTxRate       = zeros(totalRuns, T + 1);
srdkfTxRate      = zeros(totalRuns, T + 1);
rdkflocTxRate    = zeros(totalRuns, T + 1);
srdkflocTxRate   = zeros(totalRuns, T + 1);

crkfTheta            = zeros(totalRuns, T + 1);
rdkfTheta            = zeros(totalRuns, nodeCount, T + 1);
rdkfThetaBar         = zeros(totalRuns, nodeCount, T + 1);
rdkflocTheta         = zeros(totalRuns, nodeCount, T + 1);
rdkflocThetaBar      = zeros(totalRuns, nodeCount, T + 1);
srdkfTheta           = zeros(totalRuns, nodeCount, T + 1);
srdkfThetaBar        = zeros(totalRuns, nodeCount, T + 1);
srdkflocTheta        = zeros(totalRuns, nodeCount, T + 1);
srdkflocThetaBar     = zeros(totalRuns, nodeCount, T + 1);

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

  sdkfRun = sdkf.estimate(x0_hat, P0, s.X, s.Y);
  sdkfRmse(run, :)   = sdkfRun.RMSE;
  sdkfTxRate(run, :) = sdkfRun.txRate;

  srdkfRun = srdkf.estimate(x0_hat, P0, s.X, s.Y);
  srdkfRmse(run, :)   = srdkfRun.RMSE;
  srdkfTxRate(run, :) = srdkfRun.txRate;
  srdkfTheta(run, :, :) = srdkfRun.theta_hist;
  srdkfThetaBar(run, :, :) = srdkfRun.theta_bar_hist;

  rdkflocRun = rdkfloc.estimate(x0_hat, P0, s.X, s.Y);
  rdkflocRmse(run, :)   = rdkflocRun.RMSE;
  rdkflocTxRate(run, :) = rdkflocRun.txRate;
  rdkflocTheta(run, :, :) = rdkflocRun.theta_hist;
  rdkflocThetaBar(run, :, :) = rdkflocRun.theta_bar_hist;

  srdkflocRun = srdkfloc.estimate(x0_hat, P0, s.X, s.Y);
  srdkflocRmse(run, :)   = srdkflocRun.RMSE;
  srdkflocTxRate(run, :) = srdkflocRun.txRate;
  srdkflocTheta(run, :, :) = srdkflocRun.theta_hist;
  srdkflocThetaBar(run, :, :) = srdkflocRun.theta_bar_hist;
end
fprintf('Elapsed: %.3f s\n', toc);

%% Save run
results = struct( ...
  'ckfRmse', ckfRmse, 'crkfRmse', crkfRmse, 'dseacpRmse', dseacpRmse, ...
  'dkfRmse', dkfRmse, 'rdkfRmse', rdkfRmse, ...
  'sdkfRmse', sdkfRmse, 'srdkfRmse', srdkfRmse, ...
  'rdkflocRmse', rdkflocRmse, 'srdkflocRmse', srdkflocRmse, ...
  'dkfTxRate', dkfTxRate, 'rdkfTxRate', rdkfTxRate, ...
  'sdkfTxRate', sdkfTxRate, 'srdkfTxRate', srdkfTxRate, ...
  'rdkflocTxRate', rdkflocTxRate, 'srdkflocTxRate', srdkflocTxRate, ...
  'crkfTheta', crkfTheta, 'rdkfTheta', rdkfTheta, 'rdkflocTheta', rdkflocTheta, ...
  'srdkfTheta', srdkfTheta, 'srdkflocTheta', srdkflocTheta, ...
  'rdkfThetaBar', rdkfThetaBar, 'rdkflocThetaBar', rdkflocThetaBar, ...
  'srdkfThetaBar', srdkfThetaBar, 'srdkflocThetaBar', srdkflocThetaBar);
extras = struct('totalRuns', totalRuns, 'bLocal', bLocal);
savedPath = saveRun(mfilename, P, extras, netGraph, results);
