%% Simulations of the 2D Single-Target Tracking Plant under the Least-Favorable Model
%
% Mirrors estimateSST2d.m, except the per-run trajectories (X, Y) are drawn from
% the least-favorable model (Levy & Nikoukhah 2013, Section V) at tolerance b,
% rather than from the nominal plant. This is the setup used by Ghion & Zorzi
% (2023) in their Monte Carlo study.

clear;
close all;
clc;
rng(42);


%% Parameters
sst2dParams;

%% Network Definition
disp("Creating Network")
netGraph = createSpatialNetwork(nodeCount, sensorCount, maxLength);

assertConnected(netGraph);

%% Model Simulation
disp("Simulating target dynamics")
plant = SingleTarget2dModel(Ts, sensorCount, outputNoiseStd, T, turnRate);

assertStabilizable(plant.A, plant.B);
assertDetectable(plant.A, plant.C);

%% System checks (stabilizability, detectability, local observability)
disp("Checking local observability of sensor nodes")
for i = 1:sensorCount
  idx = (2 * i - 1):(2 * i);
  assertLocallyObservable(plant.A, plant.C(idx, :), i);
end
plotSystemChecks(plant.A, plant.B, plant.C, sensorCount, netGraph);

%% Estimators
disp("Initializing estimators")

% Least-favorable data generator. Forward/backward sweeps (G_t, H_t, L_t) are
% computed once here; only the per-run randn draws happen inside the parfor.
disp("Precomputing least-favorable model")
lfm = LeastFavorableModel(plant, P0, lfmKlTolerance, T);

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
%% Monte Carlo simulations
disp("Running Monte Carlo simulations (LFM data)")

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

dkfTxRate        = zeros(totalRuns, T + 1);
rdkfTxRate       = zeros(totalRuns, T + 1);
sdkfClosedTxRate = zeros(totalRuns, T + 1);
sdkfOpenTxRate   = zeros(totalRuns, T + 1);
srdkfClosedTxRate = zeros(totalRuns, T + 1);
srdkfOpenTxRate   = zeros(totalRuns, T + 1);

tic

% Showcase run (serial) - keeps sample objects available for trajectory plots.
% Also draw a nominal-model sample for side-by-side visual comparison.
mdlSample        = lfm.simulate(x0);
nomSample        = plant.simulate(x0);
ckfSample        = ckf.estimate(x0_hat, P0, mdlSample.X, mdlSample.Y);
crkfSample       = crkf.estimate(x0_hat, P0, mdlSample.X, mdlSample.Y);
dseacpSample     = dseacp.estimate(x0_hat, P0, mdlSample.X, mdlSample.Y);
dkfSample        = dkf.estimate(x0_hat, P0, mdlSample.X, mdlSample.Y);
rdkfSample       = rdkf.estimate(x0_hat, P0, mdlSample.X, mdlSample.Y);
sdkfClosedSample = sdkfClosed.estimate(x0_hat, P0, mdlSample.X, mdlSample.Y);
sdkfOpenSample   = sdkfOpen.estimate(x0_hat, P0, mdlSample.X, mdlSample.Y);
srdkfClosedSample = srdkfClosed.estimate(x0_hat, P0, mdlSample.X, mdlSample.Y);
srdkfOpenSample   = srdkfOpen.estimate(x0_hat, P0, mdlSample.X, mdlSample.Y);

ckfRmse(1, :)        = ckfSample.RMSE;
crkfRmse(1, :)       = crkfSample.RMSE;
dseacpRmse(1, :)     = dseacpSample.RMSE;
dkfRmse(1, :)        = dkfSample.RMSE;
dkfTxRate(1, :)      = dkfSample.txRate;
rdkfRmse(1, :)       = rdkfSample.RMSE;
rdkfTxRate(1, :)     = rdkfSample.txRate;
sdkfClosedRmse(1, :)   = sdkfClosedSample.RMSE;
sdkfClosedTxRate(1, :) = sdkfClosedSample.txRate;
sdkfOpenRmse(1, :)   = sdkfOpenSample.RMSE;
sdkfOpenTxRate(1, :) = sdkfOpenSample.txRate;
srdkfClosedRmse(1, :)   = srdkfClosedSample.RMSE;
srdkfClosedTxRate(1, :) = srdkfClosedSample.txRate;
srdkfOpenRmse(1, :)   = srdkfOpenSample.RMSE;
srdkfOpenTxRate(1, :) = srdkfOpenSample.txRate;

parfor run = 2:totalRuns
  fprintf('Running simulation %d/%d\n', run, totalRuns);
  s = lfm.simulate(x0);

  ckfRmse(run, :)    = ckf.estimate(x0_hat, P0, s.X, s.Y).RMSE;
  crkfRmse(run, :)   = crkf.estimate(x0_hat, P0, s.X, s.Y).RMSE;
  dseacpRmse(run, :) = dseacp.estimate(x0_hat, P0, s.X, s.Y).RMSE;

  dkfRun = dkf.estimate(x0_hat, P0, s.X, s.Y);
  dkfRmse(run, :)   = dkfRun.RMSE;
  dkfTxRate(run, :) = dkfRun.txRate;

  rdkfRun = rdkf.estimate(x0_hat, P0, s.X, s.Y);
  rdkfRmse(run, :)   = rdkfRun.RMSE;
  rdkfTxRate(run, :) = rdkfRun.txRate;

  sdkfClosedRun = sdkfClosed.estimate(x0_hat, P0, s.X, s.Y);
  sdkfClosedRmse(run, :)   = sdkfClosedRun.RMSE;
  sdkfClosedTxRate(run, :) = sdkfClosedRun.txRate;

  sdkfOpenRun = sdkfOpen.estimate(x0_hat, P0, s.X, s.Y);
  sdkfOpenRmse(run, :)   = sdkfOpenRun.RMSE;
  sdkfOpenTxRate(run, :) = sdkfOpenRun.txRate;

  srdkfClosedRun = srdkfClosed.estimate(x0_hat, P0, s.X, s.Y);
  srdkfClosedRmse(run, :)   = srdkfClosedRun.RMSE;
  srdkfClosedTxRate(run, :) = srdkfClosedRun.txRate;

  srdkfOpenRun = srdkfOpen.estimate(x0_hat, P0, s.X, s.Y);
  srdkfOpenRmse(run, :)   = srdkfOpenRun.RMSE;
  srdkfOpenTxRate(run, :) = srdkfOpenRun.txRate;
end
fprintf('Elapsed: %.3f s\n', toc);

%% Save run
results = struct( ...
  'ckfRmse', ckfRmse, 'crkfRmse', crkfRmse, 'dseacpRmse', dseacpRmse, ...
  'dkfRmse', dkfRmse, 'rdkfRmse', rdkfRmse, ...
  'sdkfClosedRmse', sdkfClosedRmse, 'sdkfOpenRmse', sdkfOpenRmse, ...
  'srdkfClosedRmse', srdkfClosedRmse, 'srdkfOpenRmse', srdkfOpenRmse, ...
  'dkfTxRate', dkfTxRate, 'rdkfTxRate', rdkfTxRate, ...
  'sdkfClosedTxRate', sdkfClosedTxRate, 'sdkfOpenTxRate', sdkfOpenTxRate, ...
  'srdkfClosedTxRate', srdkfClosedTxRate, 'srdkfOpenTxRate', srdkfOpenTxRate);
samples = struct( ...
  'mdlSample', mdlSample, 'nomSample', nomSample, ...
  'ckfSample', ckfSample, 'crkfSample', crkfSample, 'dseacpSample', dseacpSample, ...
  'dkfSample', dkfSample, 'rdkfSample', rdkfSample, ...
  'sdkfClosedSample', sdkfClosedSample, 'sdkfOpenSample', sdkfOpenSample, ...
  'srdkfClosedSample', srdkfClosedSample, 'srdkfOpenSample', srdkfOpenSample);
extras = struct('totalRuns', totalRuns);
savedPath = saveRun(mfilename, collectParams(), extras, netGraph, results, samples);

%% Plotting
plottingEnabled = true;
if plottingEnabled
  disp("Plotting results.")
  plotSST2dRun(savedPath);
end
