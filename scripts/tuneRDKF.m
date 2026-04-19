%% RDKF Hyperparameter Tuning
%
% Sweeps over the RDKF event-trigger parameters (alpha, beta, delta) plus
% the KL-divergence tolerance (b) and compares the resulting average RMSE
% and transmission rate. Based on the triggering condition from Ghion &
% Zorzi (2023) (see papers/rdkf.pdf), where a node skips its broadcast
% whenever
%
%   e_i' Omega_i e_i <= alpha   AND
%   (1/(1+beta)) Omega_i <= Psi_bar_i <= (1+delta) Omega_i  (Loewner order).
%
% A larger alpha or wider [1/(1+beta), 1+delta] range suppresses
% transmissions (saving bandwidth) at the cost of estimation accuracy. The
% KL tolerance b controls the risk-sensitivity parameter theta in the
% robust prediction step.

clear;
close all;
clc;
rng(42);

%% Fixed parameters
sst2dParams;

% totalRuns = 100;
totalRuns = 10;

%% Hyperparameter grid
alphaGrid = [0.1, 1, 5, 10];
betaGrid  = [0.1, 0.5, 1, 2.5];
deltaGrid = [0.1, 0.5, 1, 2.5];
bGrid     = [0.05, 0.1, 0.5, 1];

% Baseline values held fixed when sweeping one parameter at a time
alphaBase = 1;
betaBase  = 1;
deltaBase = 1;
bBase     = 0.05;

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
evalConfig = @(a, be, d, b) runRDKFConfig(plant, Ts, T, netGraph, ...
                                          a, be, d, b, samples, x0_hat, P0);

%% Build sweep list (one parameter varies, others held at baseline)
nConfigs = numel(alphaGrid) + numel(betaGrid) + numel(deltaGrid) + numel(bGrid);
configs = cell(nConfigs, 1);
k = 0;
for a = alphaGrid
  k = k + 1;
  configs{k} = struct('alpha', a, 'beta', betaBase, 'delta', deltaBase, ...
                      'b', bBase, 'sweep', 'alpha');
end
for be = betaGrid
  k = k + 1;
  configs{k} = struct('alpha', alphaBase, 'beta', be, 'delta', deltaBase, ...
                      'b', bBase, 'sweep', 'beta');
end
for d = deltaGrid
  k = k + 1;
  configs{k} = struct('alpha', alphaBase, 'beta', betaBase, 'delta', d, ...
                      'b', bBase, 'sweep', 'delta');
end
for bv = bGrid
  k = k + 1;
  configs{k} = struct('alpha', alphaBase, 'beta', betaBase, 'delta', deltaBase, ...
                      'b', bv, 'sweep', 'b');
end
meanRmse = zeros(nConfigs, 1);
finalRmse = zeros(nConfigs, 1);
meanTxRate = zeros(nConfigs, 1);
rmseCurves = zeros(nConfigs, T + 1);
txCurves   = zeros(nConfigs, T + 1);

disp("Running simulations...")
q = parforWaitbar(nConfigs, 'Sweeping RDKF hyperparameters');
parfor i = 1:nConfigs
  c = configs{i};
  [rmseRow, txRow] = evalConfig(c.alpha, c.beta, c.delta, c.b);
  rmseCurves(i, :) = rmseRow;
  txCurves(i, :)   = txRow;
  meanRmse(i)   = mean(rmseRow);
  finalRmse(i)  = rmseRow(end);
  meanTxRate(i) = mean(txRow);
  send(q, i);
end

%% Results table
alphaCol = arrayfun(@(k) configs{k}.alpha, 1:nConfigs)';
betaCol  = arrayfun(@(k) configs{k}.beta,  1:nConfigs)';
deltaCol = arrayfun(@(k) configs{k}.delta, 1:nConfigs)';
bCol     = arrayfun(@(k) configs{k}.b,     1:nConfigs)';
sweepCol = arrayfun(@(k) string(configs{k}.sweep), 1:nConfigs)';

resultsTable = table(sweepCol, alphaCol, betaCol, deltaCol, bCol, ...
                     meanRmse, finalRmse, meanTxRate, ...
                     'VariableNames', {'Sweep', 'alpha', 'beta', 'delta', 'b', ...
                                       'MeanRMSE', 'FinalRMSE', 'MeanTxRate'});
disp(resultsTable);

%% Plotting
disp("Plotting results")
t = (0:T) * Ts;

plotSweep(t, configs, rmseCurves, txCurves, 'alpha', ...
          {'beta', betaBase, 'delta', deltaBase, 'b', bBase});
plotSweep(t, configs, rmseCurves, txCurves, 'beta', ...
          {'alpha', alphaBase, 'delta', deltaBase, 'b', bBase});
plotSweep(t, configs, rmseCurves, txCurves, 'delta', ...
          {'alpha', alphaBase, 'beta', betaBase, 'b', bBase});
plotSweep(t, configs, rmseCurves, txCurves, 'b', ...
          {'alpha', alphaBase, 'beta', betaBase, 'delta', deltaBase});

% RMSE vs TX rate tradeoff scatter
figure;
colors = struct('alpha', [0.85 0.33 0.10], ...
                'beta',  [0.00 0.45 0.74], ...
                'delta', [0.47 0.67 0.19], ...
                'b',     [0.49 0.18 0.56]);
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
title('RDKF: RMSE vs Transmission Rate Tradeoff');
grid on;

%% Helpers

function [avgRmse, avgTx] = runRDKFConfig(plant, Ts, T, netGraph, ...
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
end

function plotSweep(t, configs, rmseCurves, txCurves, param, fixedPairs)
  idx = find(arrayfun(@(k) strcmp(configs{k}.sweep, param), 1:numel(configs)));

  fixedStr = strjoin(arrayfun(@(j) sprintf('%s=%.2g', fixedPairs{2*j-1}, fixedPairs{2*j}), ...
                              1:numel(fixedPairs)/2, 'UniformOutput', false), ', ');

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
  title(sprintf('RDKF RMSE sweep over %s (others fixed: %s)', param, fixedStr));
  legend(); grid on;

  subplot(2, 1, 2);
  hold on;
  for i = idx
    v = configs{i}.(param);
    plot(t, txCurves(i, :), 'DisplayName', sprintf('%s=%.2g', param, v));
  end
  hold off;
  xlabel('Time (s)'); ylabel('TX Rate');
  title(sprintf('RDKF Transmission Rate sweep over %s', param));
  legend(); grid on;
end
