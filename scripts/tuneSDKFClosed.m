%% SDKF (closed-loop trigger) Hyperparameter Tuning
%
% Sweeps the stochastic-trigger weight scale z, where Z = z * eye(m) and m
% is the measurement dimension. Smaller z lowers the transmission
% probability P(tx) = 1 - exp(-1/2 z'Zz), tracing an RMSE-vs-TX-rate curve.
% All other hyperparameters (alpha, beta, delta) come from sst2dParams.m; in
% an all-sensor network they are unused, since every node runs the
% stochastic trigger rather than the deterministic relay trigger.
%
% Saved via saveRun and plotted by plotTuneRun. For a matched-rate
% comparison against DKF, run this and tuneDKF, then overlay:
%   plotTuneRun('results/tuneSDKFClosed_...mat', 'results/tuneDKF_...mat')
% and set sensorCount = nodeCount in sst2dParams so the stochastic trigger
% governs every node (otherwise the deterministic relay trigger dominates).

clear;
close all;
clc;
rng(42);

%% Fixed parameters
sst2dParams;

%% Hyperparameter grid
zGrid = [1e-5, 2e-5, 3e-5, 4e-5, 4.5e-5, 5e-5, 6e-5, 7e-5, 1e-4];

%% Network and Plant
disp("Creating Network")
netGraph = createSpatialNetwork(nodeCount, sensorCount, maxLength);

disp("Simulating target dynamics")
plant = SingleTarget2dModel(Ts, sensorCount, outputNoiseStd, T, turnRate);

%% Pre-generate trajectories so all configurations see the same data
disp("Pre-generating Monte Carlo trajectories")
samples = cell(totalRuns, 1);
for run = 1:totalRuns
  samples{run} = plant.simulate(x0);
end

%% Sweep
zDim = 2; % measurement dimension (m), Z in S^m_{++} per Han et al. 2015
nConfigs = numel(zGrid);
configs = cell(nConfigs, 1);
for k = 1:nConfigs
  configs{k} = struct('z', zGrid(k), 'sweep', 'z');
end

meanRmse   = zeros(nConfigs, 1);
finalRmse  = zeros(nConfigs, 1);
meanTxRate = zeros(nConfigs, 1);
ssRmseMean = zeros(nConfigs, 1);
ssRmseStd  = zeros(nConfigs, 1);
rmseCurves = zeros(nConfigs, T + 1);
txCurves   = zeros(nConfigs, T + 1);

disp("Running simulations...")
parfor i = 1:nConfigs
  fprintf('Sweeping SDKF-Closed Z %d/%d\n', i, nConfigs);
  c = configs{i};
  [rmseRow, txRow, ssM, ssS] = runSDKFConfig(plant, Ts, T, netGraph, ...
                                             dkfAlpha, dkfBeta, dkfDelta, ...
                                             c.z * eye(zDim), 'closed', ...
                                             samples, x0_hat, P0);
  rmseCurves(i, :) = rmseRow;
  txCurves(i, :)   = txRow;
  meanRmse(i)   = mean(rmseRow);
  finalRmse(i)  = rmseRow(end);
  meanTxRate(i) = mean(txRow(2:end));   % drop the t=0 sample (always 0)
  ssRmseMean(i) = ssM;
  ssRmseStd(i)  = ssS;
end

%% Results table
zCol     = arrayfun(@(k) configs{k}.z, 1:nConfigs)';
sweepCol = arrayfun(@(k) string(configs{k}.sweep), 1:nConfigs)';

resultsTable = table(sweepCol, zCol, meanRmse, finalRmse, meanTxRate, ssRmseMean, ssRmseStd, ...
                     'VariableNames', {'Sweep', 'z', ...
                                       'MeanRMSE', 'FinalRMSE', 'MeanTxRate', ...
                                       'SsRMSE', 'SsRMSEStd'});
disp(resultsTable);

%% Save run
results = struct( ...
  'rmseCurves', rmseCurves, 'txCurves', txCurves, ...
  'meanRmse', meanRmse, 'finalRmse', finalRmse, 'meanTxRate', meanTxRate, ...
  'ssRmseMean', ssRmseMean, 'ssRmseStd', ssRmseStd, ...
  'configs', {configs}, 'resultsTable', resultsTable);
extras = struct( ...
  'totalRuns', totalRuns, 'filterName', 'SDKF-Closed', ...
  'zGrid', zGrid, ...
  'bases', struct('z', NaN));
savedPath = saveRun(mfilename, collectParams(), extras, netGraph, results, struct());

%% Plotting
disp("Plotting results")
plotTuneRun(savedPath);

%% Helpers

function [avgRmse, avgTx, ssMean, ssStd] = runSDKFConfig(plant, Ts, T, netGraph, ...
                                                         alpha, beta, delta, Z, triggerMode, ...
                                                         samples, x0_hat, P0)
  sdkf = SDKF(plant, Ts, T, netGraph, alpha, beta, delta, Z, triggerMode);
  nRuns = numel(samples);
  rmseLog = zeros(nRuns, T + 1);
  txLog   = zeros(nRuns, T + 1);
  for run = 1:nRuns
    s = samples{run};
    out = sdkf.estimate(x0_hat, P0, s.X, s.Y);
    rmseLog(run, :) = out.RMSE;
    txLog(run, :)   = out.txRate;
  end
  avgRmse = mean(rmseLog, 1);
  avgTx   = mean(txLog, 1);
  [ssMean, ssStd] = ssRmseStats(rmseLog, T);
end
