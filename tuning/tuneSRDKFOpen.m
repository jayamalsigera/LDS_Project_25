%% SRDKF (open-loop trigger) Hyperparameter Tuning
%
% Sweeps the error-norm weight scale z, where Z = z * eye(2). The
% open-loop trigger uses only the measured output (hence Z is 2 x 2 for
% the 2D target). All other hyperparameters (alpha, beta, delta, KL
% tolerance) come from sst2dParams.m.

clear;
close all;
clc;
rng(42);

%% Fixed parameters
sst2dParams;

totalTuneRuns = 50;   % fewer runs than the estimate scripts; enough to separate configs

%% Hyperparameter grid
zGrid = [0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300];

%% Network and Plant
disp("Creating Network")
netGraph = createSpatialNetwork(nodeCount, sensorCount, maxLength);

disp("Simulating target dynamics")
plant = SingleTarget2dModel(Ts, sensorCount, outputNoiseStd, T, turnRate);

%% Pre-generate trajectories so all configurations see the same data
disp("Pre-generating Monte Carlo trajectories")
samples = cell(totalTuneRuns, 1);
for run = 1:totalTuneRuns
  samples{run} = plant.simulate(x0);
end

%% Sweep
zDim = 2;
nConfigs = numel(zGrid);
configs = cell(nConfigs, 1);
for k = 1:nConfigs
  configs{k} = struct('z', zGrid(k), 'sweep', 'z');
end

disp("Running simulations...")
makeFilter = @(c) SRDKF(plant, Ts, T, netGraph, dkfAlpha, dkfBeta, dkfDelta, ...
                        klTolerance, c.z * eye(zDim), 'open');
[meanRmse, finalRmse, meanTxRate, ssRmseMean, ssRmseStd, rmseCurves, txCurves] = ...
    evalConfigsMC(makeFilter, configs, samples, x0_hat, P0, T);

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
  'totalRuns', totalTuneRuns, 'filterName', 'SRDKF-Open', ...
  'zGrid', zGrid, ...
  'bases', struct('z', NaN));
savedPath = saveRun(mfilename, collectParams(), extras, netGraph, results, struct());

%% Plotting
disp("Plotting results")
plotTuneRun(savedPath);
