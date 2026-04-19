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

%% Estimators
disp("Initializing estimators")

% Least-favorable data generator. Forward/backward sweeps (G_t, H_t, L_t) are
% computed once here; only the per-run randn draws happen inside the parfor.
disp("Precomputing least-favorable model")
lfm = LeastFavorableModel(plant, P0, klTolerance, T);

% Initiate all filters
ckf = CKF(plant, Ts, T);
dseacp = DSEACP(plant, Ts, T, netGraph, consensusSteps);
dkf = DKF(plant, Ts, T, netGraph, dkfAlpha, dkfBeta, dkfDelta);
rdkf = RDKF(plant, Ts, T, netGraph, dkfAlpha, dkfBeta, dkfDelta, klTolerance);
srdkfClosed = SRDKF(plant, Ts, T, netGraph, dkfAlpha, dkfBeta, dkfDelta, ...
                    klTolerance, errorNormWeightsClosed, 'closed');
srdkfOpen = SRDKF(plant, Ts, T, netGraph, dkfAlpha, dkfBeta, dkfDelta, ...
                  klTolerance, errorNormWeightsOpen, 'open');
%% Monte Carlo simulations
disp("Running Monte Carlo simulations (LFM data)")

% totalRuns = 200;
totalRuns = 100;
% totalRuns = 10;
% totalRuns = 2;

% Preallocate RMSE and Transmission rate logs
ckfRmse = zeros(totalRuns, T + 1);
dseacpRmse = zeros(totalRuns, T + 1);
dkfRmse = zeros(totalRuns, T + 1);
rdkfRmse = zeros(totalRuns, T + 1);
srdkfClosedRmse = zeros(totalRuns, T + 1);
srdkfOpenRmse   = zeros(totalRuns, T + 1);

dkfTxRate = zeros(totalRuns, T + 1);
rdkfTxRate = zeros(totalRuns, T + 1);
srdkfClosedTxRate = zeros(totalRuns, T + 1);
srdkfOpenTxRate   = zeros(totalRuns, T + 1);

tic

% Showcase run (serial) — keeps sample objects available for trajectory plots.
% Also draw a nominal-model sample for side-by-side visual comparison.
mdlSample = lfm.simulate(x0);
nomSample = plant.simulate(x0);
ckfSample = ckf.estimate(x0_hat, P0, mdlSample.X, mdlSample.Y);
dseacpSample = dseacp.estimate(x0_hat, P0, mdlSample.X, mdlSample.Y);
dkfSample = dkf.estimate(x0_hat, P0, mdlSample.X, mdlSample.Y);
rdkfSample = rdkf.estimate(x0_hat, P0, mdlSample.X, mdlSample.Y);
srdkfClosedSample = srdkfClosed.estimate(x0_hat, P0, mdlSample.X, mdlSample.Y);
srdkfOpenSample = srdkfOpen.estimate(x0_hat, P0, mdlSample.X, mdlSample.Y);

ckfRmse(1, :) = ckfSample.RMSE;
dseacpRmse(1, :) = dseacpSample.RMSE;
dkfRmse(1, :) = dkfSample.RMSE;
dkfTxRate(1, :) = dkfSample.txRate;
rdkfRmse(1, :) = rdkfSample.RMSE;
rdkfTxRate(1, :) = rdkfSample.txRate;
srdkfClosedRmse(1, :) = srdkfClosedSample.RMSE;
srdkfClosedTxRate(1, :) = srdkfClosedSample.txRate;
srdkfOpenRmse(1, :) = srdkfOpenSample.RMSE;
srdkfOpenTxRate(1, :) = srdkfOpenSample.txRate;

q = parforWaitbar(totalRuns - 1, 'Running simulations');
parfor run = 2:totalRuns
  % Draw one trajectory from the least-favorable model
  s = lfm.simulate(x0);

  % Run each estimator on the same data and log RMSE + transmission rate
  ckfRmse(run, :) = ckf.estimate(x0_hat, P0, s.X, s.Y).RMSE;
  dseacpRmse(run, :) = dseacp.estimate(x0_hat, P0, s.X, s.Y).RMSE;

  dkfRun = dkf.estimate(x0_hat, P0, s.X, s.Y);
  dkfRmse(run, :) = dkfRun.RMSE;
  dkfTxRate(run, :) = dkfRun.txRate;

  rdkfRun = rdkf.estimate(x0_hat, P0, s.X, s.Y);
  rdkfRmse(run, :) = rdkfRun.RMSE;
  rdkfTxRate(run, :) = rdkfRun.txRate;

  srdkfClosedRun = srdkfClosed.estimate(x0_hat, P0, s.X, s.Y);
  srdkfClosedRmse(run, :) = srdkfClosedRun.RMSE;
  srdkfClosedTxRate(run, :) = srdkfClosedRun.txRate;

  srdkfOpenRun = srdkfOpen.estimate(x0_hat, P0, s.X, s.Y);
  srdkfOpenRmse(run, :) = srdkfOpenRun.RMSE;
  srdkfOpenTxRate(run, :) = srdkfOpenRun.txRate;

  send(q, run);
end
fprintf('Elapsed: %.3f s\n', toc);

%% Save run
results = struct( ...
  'ckfRmse', ckfRmse, 'dseacpRmse', dseacpRmse, ...
  'dkfRmse', dkfRmse, 'rdkfRmse', rdkfRmse, ...
  'srdkfClosedRmse', srdkfClosedRmse, 'srdkfOpenRmse', srdkfOpenRmse, ...
  'dkfTxRate', dkfTxRate, 'rdkfTxRate', rdkfTxRate, ...
  'srdkfClosedTxRate', srdkfClosedTxRate, 'srdkfOpenTxRate', srdkfOpenTxRate);
samples = struct( ...
  'mdlSample', mdlSample, 'nomSample', nomSample, ...
  'ckfSample', ckfSample, 'dseacpSample', dseacpSample, ...
  'dkfSample', dkfSample, 'rdkfSample', rdkfSample, ...
  'srdkfClosedSample', srdkfClosedSample, 'srdkfOpenSample', srdkfOpenSample);
extras = struct('totalRuns', totalRuns);
savedPath = saveRun(mfilename, collectParams(), extras, netGraph, results, samples);

%% Plotting
plottingEnabled = true;
if plottingEnabled
  disp("Plotting results.")
  plotSST2dRun(loadRun(savedPath));
end
