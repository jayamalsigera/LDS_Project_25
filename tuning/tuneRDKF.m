%% RDKF Hyperparameter Tuning (delta x beta grid)
%
% Full-factorial sweep over the RDKF event-trigger parameters delta and beta
% (alpha and the KL tolerance b held fixed), based on the triggering
% condition from Ghion & Zorzi (2023) (see papers/rdkf.pdf), where a node
% skips its broadcast whenever
%
%   e_i' Omega_i e_i <= alpha   AND
%   (1/(1+beta)) Omega_i <= Psi_bar_i <= (1+delta) Omega_i  (Loewner order).
%
% As with DKF, one-at-a-time sweeps miss the interaction between the two
% Loewner-tolerance knobs (beta, delta), so this grid measures every corner
% and extracts the Pareto frontier of the RMSE-vs-TX tradeoff. alpha and the
% robust-prediction tolerance b are fixed; b in particular is the robust
% layer's knob, not a trigger knob, and is left at its baseline here.

clear;
close all;
clc;
rng(42);

%% Fixed parameters
sst2dParams;

totalTuneRuns = 50;   % fewer runs than the estimate scripts; enough to separate configs

%% Hyperparameter grid (delta x beta; alpha and b fixed)
betaGrid  = [0.1, 0.5, 1, 2.5];
deltaGrid = [0.1, 0.5, 1, 2.5];
alphaFixed = 1;
bFixed     = 0.05;

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

%% Build grid (beta outer, delta inner -> reshape(v, nDelta, nBeta))
nBeta  = numel(betaGrid);
nDelta = numel(deltaGrid);
nConfigs = nBeta * nDelta;
configs = cell(nConfigs, 1);
k = 0;
for bi = 1:nBeta
  for di = 1:nDelta
    k = k + 1;
    configs{k} = struct('alpha', alphaFixed, 'beta', betaGrid(bi), ...
                        'delta', deltaGrid(di), 'b', bFixed, 'sweep', 'grid');
  end
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
  fprintf('Sweeping RDKF grid %d/%d\n', i, nConfigs);
  c = configs{i};
  [rmseRow, txRow, ssM, ssS] = runRDKFConfig(plant, Ts, T, netGraph, ...
                                             c.alpha, c.beta, c.delta, c.b, ...
                                             samples, x0_hat, P0);
  rmseCurves(i, :) = rmseRow;
  txCurves(i, :)   = txRow;
  meanRmse(i)   = mean(rmseRow);
  finalRmse(i)  = rmseRow(end);
  meanTxRate(i) = mean(txRow(2:end));   % drop the t=0 sample (always 0)
  ssRmseMean(i) = ssM;
  ssRmseStd(i)  = ssS;
end

%% Results table + Pareto frontier
betaCol  = arrayfun(@(k) configs{k}.beta,  1:nConfigs)';
deltaCol = arrayfun(@(k) configs{k}.delta, 1:nConfigs)';
onFront  = paretoFront(meanTxRate, ssRmseMean);

resultsTable = table(betaCol, deltaCol, meanRmse, finalRmse, meanTxRate, ...
                     ssRmseMean, ssRmseStd, onFront, ...
                     'VariableNames', {'beta', 'delta', ...
                                       'MeanRMSE', 'FinalRMSE', 'MeanTxRate', ...
                                       'SsRMSE', 'SsRMSEStd', 'Pareto'});
resultsTable = sortrows(resultsTable, 'MeanTxRate');
disp(resultsTable);
disp('Pareto-optimal configs (min TX, min RMSE):');
disp(sortrows(resultsTable(resultsTable.Pareto, :), 'MeanTxRate'));

%% Save run
results = struct( ...
  'rmseCurves', rmseCurves, 'txCurves', txCurves, ...
  'meanRmse', meanRmse, 'finalRmse', finalRmse, 'meanTxRate', meanTxRate, ...
  'ssRmseMean', ssRmseMean, 'ssRmseStd', ssRmseStd, 'onFront', onFront, ...
  'configs', {configs}, 'resultsTable', resultsTable);
extras = struct( ...
  'totalRuns', totalTuneRuns, 'filterName', 'RDKF', ...
  'alphaFixed', alphaFixed, 'bFixed', bFixed, ...
  'gridAxisCol', struct('name', 'beta',  'values', betaGrid), ...   % reshape columns
  'gridAxisRow', struct('name', 'delta', 'values', deltaGrid));     % reshape rows
savedPath = saveRun(mfilename, collectParams(), extras, netGraph, results, struct());

%% Plotting
disp("Plotting results")
plotTuneRun(savedPath);

%% Helpers

function [avgRmse, avgTx, ssMean, ssStd] = runRDKFConfig(plant, Ts, T, netGraph, ...
                                                         alpha, beta, delta, b, samples, x0_hat, P0)
  rdkf = RDKF(plant, Ts, T, netGraph, alpha, beta, delta, b);
  nRuns = numel(samples);
  rmseLog = zeros(nRuns, T + 1);
  txLog   = zeros(nRuns, T + 1);
  for run = 1:nRuns
    s = samples{run};
    out = rdkf.estimate(x0_hat, P0, s.X, s.Y);
    rmseLog(run, :) = out.RMSE;
    txLog(run, :)   = out.txRate;
  end
  avgRmse = mean(rmseLog, 1);
  avgTx   = mean(txLog, 1);
  [ssMean, ssStd] = ssRmseStats(rmseLog, T);
end
