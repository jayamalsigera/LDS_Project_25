function [ssMean, ssStd] = ssRmseStats(rmseLog, T)
%SSRMSESTATS  Steady-state RMSE mean and std across Monte Carlo runs.
%
% The steady-state window is the second half of the horizon, so transient
% convergence does not pollute the summary. rmseLog is nRuns-by-(T+1); the
% per-run steady-state RMSE is averaged over that window, then reduced to a
% mean and standard deviation across runs (the std feeds plot error bars).

  ssWin    = (round(T / 2) + 1):(T + 1);
  ssPerRun = mean(rmseLog(:, ssWin), 2);
  ssMean   = mean(ssPerRun);
  ssStd    = std(ssPerRun);
end
