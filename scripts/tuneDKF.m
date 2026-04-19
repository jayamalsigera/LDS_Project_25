%% DKF Hyperparameter Tuning
%
% Sweeps over the DKF event-trigger parameters (alpha, beta, delta) and
% compares the resulting average RMSE and transmission rate. Based on the
% triggering condition from Battistelli et al. (see papers/dkf.pdf), where
% a node broadcasts its local innovation whenever
%
%   ||y_i - C_i x_hat_i|| > alpha * beta^k + delta.
%
% A larger alpha/delta or smaller beta suppresses transmissions (saving
% bandwidth) at the cost of estimation accuracy.

clear;
close all;
clc;
rng(42);

%% Fixed parameters
sst2dParams;

%% Hyperparameter grid
alphaGrid = [0.01, 0.1, 1, 10];
betaGrid  = [0.1, 0.2, 0.5, 0.9];
deltaGrid = [0.1, 0.5, 1.0];

% Baseline values held fixed when sweeping one parameter at a time
alphaBase = 1;
betaBase  = 0.2;
deltaBase = 0.5;

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

%% Sweep helper
evalConfig = @(a, b, d) runDKFConfig(plant, Ts, T, netGraph, ...
                                     a, b, d, samples, x0_hat, P0);

%% Build sweep list (one parameter varies, other two held at baseline)
nConfigs = numel(alphaGrid) + numel(betaGrid) + numel(deltaGrid);
configs = cell(nConfigs, 1);
k = 0;
for a = alphaGrid
  k = k + 1;
  configs{k} = struct('alpha', a, 'beta', betaBase, 'delta', deltaBase, ...
                      'sweep', 'alpha');
end
for b = betaGrid
  k = k + 1;
  configs{k} = struct('alpha', alphaBase, 'beta', b, 'delta', deltaBase, ...
                      'sweep', 'beta');
end
for d = deltaGrid
  k = k + 1;
  configs{k} = struct('alpha', alphaBase, 'beta', betaBase, 'delta', d, ...
                      'sweep', 'delta');
end
meanRmse = zeros(nConfigs, 1);
finalRmse = zeros(nConfigs, 1);
meanTxRate = zeros(nConfigs, 1);
rmseCurves = zeros(nConfigs, T + 1);
txCurves   = zeros(nConfigs, T + 1);

disp("Running simulations...")
parfor i = 1:nConfigs
  fprintf('Sweeping DKF hyperparameters %d/%d\n', i, nConfigs);
  c = configs{i};
  [rmseRow, txRow] = evalConfig(c.alpha, c.beta, c.delta);
  rmseCurves(i, :) = rmseRow;
  txCurves(i, :)   = txRow;
  meanRmse(i)   = mean(rmseRow);
  finalRmse(i)  = rmseRow(end);
  meanTxRate(i) = mean(txRow);
end

%% Results table
alphaCol = arrayfun(@(k) configs{k}.alpha, 1:nConfigs)';
betaCol  = arrayfun(@(k) configs{k}.beta,  1:nConfigs)';
deltaCol = arrayfun(@(k) configs{k}.delta, 1:nConfigs)';
sweepCol = arrayfun(@(k) string(configs{k}.sweep), 1:nConfigs)';

resultsTable = table(sweepCol, alphaCol, betaCol, deltaCol, ...
                     meanRmse, finalRmse, meanTxRate, ...
                     'VariableNames', {'Sweep', 'alpha', 'beta', 'delta', ...
                                       'MeanRMSE', 'FinalRMSE', 'MeanTxRate'});
disp(resultsTable);

%% Save run
results = struct( ...
  'rmseCurves', rmseCurves, 'txCurves', txCurves, ...
  'meanRmse', meanRmse, 'finalRmse', finalRmse, 'meanTxRate', meanTxRate, ...
  'configs', {configs}, 'resultsTable', resultsTable);
extras = struct( ...
  'totalRuns', totalRuns, ...
  'alphaGrid', alphaGrid, 'betaGrid', betaGrid, 'deltaGrid', deltaGrid, ...
  'alphaBase', alphaBase, 'betaBase', betaBase, 'deltaBase', deltaBase);
savedPath = saveRun(mfilename, collectParams(), extras, netGraph, results, struct());

%% Plotting
disp("Plotting results")
plotTuneRun(loadRun(savedPath));

%% Helpers

function [avgRmse, avgTx] = runDKFConfig(plant, Ts, T, netGraph, ...
                                         alpha, beta, delta, samples, x0_hat, P0)
  dkf = DKF(plant, Ts, T, netGraph, alpha, beta, delta);
  nRuns = numel(samples);
  rmseLog = zeros(nRuns, T + 1);
  txLog   = zeros(nRuns, T + 1);
  for run = 1:nRuns
    s = samples{run};
    out = dkf.estimate(x0_hat, P0, s.X, s.Y);
    rmseLog(run, :) = out.RMSE;
    txLog(run, :)   = out.txRate;
  end
  avgRmse = mean(rmseLog, 1);
  avgTx   = mean(txLog, 1);
end
