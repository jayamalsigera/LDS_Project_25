%% SRDKF (closed-loop trigger) Hyperparameter Tuning
%
% Sweeps the error-norm weight scale z, where Z = z * eye(length(x0)).
% All other hyperparameters (alpha, beta, delta, KL tolerance) come from
% sst2dParams.m.

clear;
close all;
clc;
rng(42);

%% Fixed parameters
sst2dParams;

%% Hyperparameter grid
zGrid = [0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300];

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
zDim = length(x0);
nConfigs = numel(zGrid);
configs = cell(nConfigs, 1);
for k = 1:nConfigs
  configs{k} = struct('z', zGrid(k), 'sweep', 'z');
end

meanRmse = zeros(nConfigs, 1);
finalRmse = zeros(nConfigs, 1);
meanTxRate = zeros(nConfigs, 1);
rmseCurves = zeros(nConfigs, T + 1);
txCurves   = zeros(nConfigs, T + 1);

disp("Running simulations...")
parfor i = 1:nConfigs
  fprintf('Sweeping SRDKF-Closed Z %d/%d\n', i, nConfigs);
  c = configs{i};
  [rmseRow, txRow] = runSRDKFConfig(plant, Ts, T, netGraph, ...
                                    dkfAlpha, dkfBeta, dkfDelta, klTolerance, ...
                                    c.z * eye(zDim), 'closed', ...
                                    samples, x0_hat, P0);
  rmseCurves(i, :) = rmseRow;
  txCurves(i, :)   = txRow;
  meanRmse(i)   = mean(rmseRow);
  finalRmse(i)  = rmseRow(end);
  meanTxRate(i) = mean(txRow);
end

%% Results table
zCol     = arrayfun(@(k) configs{k}.z, 1:nConfigs)';
sweepCol = arrayfun(@(k) string(configs{k}.sweep), 1:nConfigs)';

resultsTable = table(sweepCol, zCol, meanRmse, finalRmse, meanTxRate, ...
                     'VariableNames', {'Sweep', 'z', ...
                                       'MeanRMSE', 'FinalRMSE', 'MeanTxRate'});
disp(resultsTable);

%% Save run
results = struct( ...
  'rmseCurves', rmseCurves, 'txCurves', txCurves, ...
  'meanRmse', meanRmse, 'finalRmse', finalRmse, 'meanTxRate', meanTxRate, ...
  'configs', {configs}, 'resultsTable', resultsTable);
extras = struct( ...
  'totalRuns', totalRuns, 'filterName', 'SRDKF-Closed', ...
  'zGrid', zGrid, ...
  'bases', struct('z', NaN));
savedPath = saveRun(mfilename, collectParams(), extras, netGraph, results, struct());

%% Plotting
disp("Plotting results")
plotTuneRun(savedPath);

%% Helpers

function [avgRmse, avgTx] = runSRDKFConfig(plant, Ts, T, netGraph, ...
                                           alpha, beta, delta, b, Z, triggerMode, ...
                                           samples, x0_hat, P0)
  srdkf = SRDKF(plant, Ts, T, netGraph, alpha, beta, delta, b, Z, triggerMode);
  nRuns = numel(samples);
  rmseLog = zeros(nRuns, T + 1);
  txLog   = zeros(nRuns, T + 1);
  for run = 1:nRuns
    s = samples{run};
    out = srdkf.estimate(x0_hat, P0, s.X, s.Y);
    rmseLog(run, :) = out.RMSE;
    txLog(run, :)   = out.txRate;
  end
  avgRmse = mean(rmseLog, 1);
  avgTx   = mean(txLog, 1);
end
