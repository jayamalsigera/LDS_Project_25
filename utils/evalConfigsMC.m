function [meanRmse, finalRmse, meanTx, ssMean, ssStd, rmseCurves, txCurves] = ...
    evalConfigsMC(makeFilter, configs, samples, x0_hat, P0, T)
%EVALCONFIGSMC  Monte Carlo evaluation of a set of estimator configs.
%
%   [...] = evalConfigsMC(makeFilter, configs, samples, x0_hat, P0, T)
%   evaluates every config in CONFIGS against every trajectory in SAMPLES.
%
%   The work is flattened into a single parfor over all (config, run) pairs
%   rather than a parfor over configs with the runs looping serially inside.
%   This raises parallel width from numel(configs) to numel(configs) *
%   numel(samples), so cluster CPU efficiency is no longer capped at
%   nConfigs/nCores, and the fine granularity averages out per-config runtime
%   imbalance. Common random numbers are preserved: run ri always uses
%   samples{ri}, so every config sees identical trajectories.
%
%   MAKEFILTER is a handle @(cfg) -> estimator object; it closes over the
%   plant, Ts, T and network and builds a fresh filter per (config, run) pair
%   inside the worker (so no estimator object is broadcast).
%
%   Returns per-config column vectors (meanRmse, finalRmse, meanTx, ssMean,
%   ssStd) and per-config averaged curves (rmseCurves, txCurves), matching the
%   arrays the tune* scripts previously built by hand.

  setupPool();   % match the pool to the SLURM allocation (no-op off-cluster)

  nC = numel(configs);
  nR = numel(samples);
  P  = nC * nR;

  allRmse = zeros(P, T + 1);
  allTx   = zeros(P, T + 1);

  parfor p = 1:P
    ci = floor((p - 1) / nR) + 1;
    ri = mod(p - 1, nR) + 1;
    f  = makeFilter(configs{ci});          % built inside the worker
    s  = samples{ri};
    out = f.estimate(x0_hat, P0, s.X, s.Y);
    allRmse(p, :) = out.RMSE;
    allTx(p, :)   = out.txRate;
  end

  % Reduce: rows [ (ci-1)*nR + (1:nR) ] belong to config ci.
  meanRmse   = zeros(nC, 1); finalRmse = zeros(nC, 1); meanTx = zeros(nC, 1);
  ssMean     = zeros(nC, 1); ssStd     = zeros(nC, 1);
  rmseCurves = zeros(nC, T + 1); txCurves = zeros(nC, T + 1);
  for ci = 1:nC
    rows = (ci - 1) * nR + (1:nR);
    rLog = allRmse(rows, :);
    tLog = allTx(rows, :);
    rmseCurves(ci, :) = mean(rLog, 1);
    txCurves(ci, :)   = mean(tLog, 1);
    meanRmse(ci)  = mean(rmseCurves(ci, :));
    finalRmse(ci) = rmseCurves(ci, end);
    meanTx(ci)    = mean(txCurves(ci, 2:end));   % drop t=0 (always 0)
    [ssMean(ci), ssStd(ci)] = ssRmseStats(rLog, T);
  end
end
