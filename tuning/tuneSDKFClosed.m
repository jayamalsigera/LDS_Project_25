function tuneSDKFClosed(scenario)
%% SDKF (closed-loop trigger) Hyperparameter Tuning
%
%   tuneSDKFClosed('sst2d')   tuneSDKFClosed('sst3d')
%
% Sweeps the stochastic-trigger weight scale z, where Z = z * eye(m) and m
% is the measurement dimension. Smaller z lowers the transmission
% probability P(tx) = 1 - exp(-1/2 z'Zz), tracing an RMSE-vs-TX-rate curve.
% All other hyperparameters (alpha, beta, delta) come from the scenario's
% params file; in an all-sensor network they are unused, since every node
% runs the stochastic trigger rather than the deterministic relay trigger.
%
% Saved via saveRun and plotted by plotTuneRun. For a matched-rate
% comparison against DKF, run this and tuneDKF, then overlay:
%   plotTuneRun('results/tuneSDKFClosed_...mat', 'results/tuneDKF_...mat')
% (keep sensorCount = nodeCount so the stochastic trigger governs every node,
% otherwise the deterministic relay trigger dominates).

  rng(42);

  %% Fixed parameters
  P = sstTuneParams(scenario);
  Ts = P.Ts; T = P.T;
  totalTuneRuns = P.totalTuneRuns;

  %% Hyperparameter grid
  zGrid = P.tuneSdkfZGrid;

  %% Network and Plant
  disp("Creating Network")
  netGraph = createSpatialNetwork(P.nodeCount, P.sensorCount, P.maxLength);

  disp("Simulating target dynamics")
  [plant, zDim] = sstPlant(P);   % zDim = measurement dimension m

  %% Pre-generate trajectories so all configurations see the same data
  disp("Pre-generating Monte Carlo trajectories")
  samples = cell(totalTuneRuns, 1);
  for run = 1:totalTuneRuns
    samples{run} = plant.simulate(P.x0);
  end

  %% Sweep
  nConfigs = numel(zGrid);
  configs = cell(nConfigs, 1);
  for k = 1:nConfigs
    configs{k} = struct('z', zGrid(k), 'sweep', 'z');
  end

  disp("Running simulations...")
  makeFilter = @(c) SDKF(plant, Ts, T, netGraph, P.dkfDelta, c.z * eye(zDim), 'closed');
  [meanRmse, finalRmse, meanTxRate, ssRmseMean, ssRmseStd, rmseCurves, txCurves] = ...
      evalConfigsMC(makeFilter, configs, samples, P.x0_hat, P.P0, T);

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
    'totalRuns', totalTuneRuns, 'filterName', 'SDKF-Closed', ...
    'zGrid', zGrid, ...
    'bases', struct('z', NaN));
  savedPath = saveRun(sprintf('tuneSDKFClosed_%s', scenario), P, extras, netGraph, results, struct());

  %% Plotting
  disp("Plotting results")
  plotTuneRun(savedPath);
end
