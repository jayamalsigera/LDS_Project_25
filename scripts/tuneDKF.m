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
T = 1000;
Ts = 0.1;
outputNoiseStd = 10;

x0 = [25 25 1000 1000]';
nodeCount = 20;
sensorCount = 4;
maxLength = 5000;

totalRuns = 100;
% totalRuns = 5;

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
plant = SingleTarget2dModel(Ts, sensorCount, outputNoiseStd, T);

x0_hat = x0;
P0 = 1e3 * eye(size(x0, 1));

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
  c = configs{i};
  tConfig = tic;
  [rmseRow, txRow] = evalConfig(c.alpha, c.beta, c.delta);
  rmseCurves(i, :) = rmseRow;
  txCurves(i, :)   = txRow;
  meanRmse(i)   = mean(rmseRow);
  finalRmse(i)  = rmseRow(end);
  meanTxRate(i) = mean(txRow);
  fprintf('Config %d/%d done (%.3f s)\n', i, nConfigs, toc(tConfig));
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

%% Plotting
disp("Plotting results")
t = (0:T) * Ts;

plotSweep(t, configs, rmseCurves, txCurves, 'alpha', betaBase, deltaBase);
plotSweep(t, configs, rmseCurves, txCurves, 'beta',  alphaBase, deltaBase);
plotSweep(t, configs, rmseCurves, txCurves, 'delta', alphaBase, betaBase);

% RMSE vs TX rate tradeoff scatter
figure;
colors = struct('alpha', [0.85 0.33 0.10], ...
                'beta',  [0.00 0.45 0.74], ...
                'delta', [0.47 0.67 0.19]);
hold on;
for i = 1:nConfigs
  c = configs{i};
  col = colors.(c.sweep);
  scatter(meanTxRate(i), meanRmse(i), 60, col, 'filled', ...
          'DisplayName', sprintf('%s=%.2g', c.sweep, c.(c.sweep)));
  text(meanTxRate(i), meanRmse(i), ...
       sprintf(' %s=%.2g', c.sweep, c.(c.sweep)), 'FontSize', 8);
end
hold off;
xlabel('Mean Transmission Rate');
ylabel('Mean RMSE');
title('DKF: RMSE vs Transmission Rate Tradeoff');
grid on;

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

function plotSweep(t, configs, rmseCurves, txCurves, param, fixedA, fixedB)
  idx = find(arrayfun(@(k) strcmp(configs{k}.sweep, param), 1:numel(configs)));

  figure;
  subplot(2, 1, 1);
  hold on;
  for i = idx
    v = configs{i}.(param);
    semilogy(t, rmseCurves(i, :), 'DisplayName', sprintf('%s=%.2g', param, v));
  end
  hold off;
  set(gca, 'YScale', 'log');
  xlabel('Time (s)'); ylabel('RMSE');
  title(sprintf('DKF RMSE sweep over %s (others fixed: %.2g, %.2g)', ...
                param, fixedA, fixedB));
  legend(); grid on;

  subplot(2, 1, 2);
  hold on;
  for i = idx
    v = configs{i}.(param);
    plot(t, txCurves(i, :), 'DisplayName', sprintf('%s=%.2g', param, v));
  end
  hold off;
  xlabel('Time (s)'); ylabel('TX Rate');
  title(sprintf('DKF Transmission Rate sweep over %s', param));
  legend(); grid on;
end
