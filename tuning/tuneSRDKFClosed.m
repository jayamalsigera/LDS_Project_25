function tuneSRDKFClosed(scenario)
%% SRDKF (closed-loop trigger) Hyperparameter Tuning
%
%   tuneSRDKFClosed('sst2d')   tuneSRDKFClosed('sst3d')
%
% Sweeps the error-norm weight scale z, where Z = z * eye(m) and m is the
% measurement dimension. All other hyperparameters (alpha, beta, delta, KL
% tolerance) come from the scenario's params file. The z grid lives there
% (tuneSrdkfZGrid).

  rng(42);

  %% Fixed parameters
  P = sstTuneParams(scenario);
  Ts = P.Ts; T = P.T;
  totalTuneRuns = P.totalTuneRuns;

  %% Hyperparameter grid
  zGrid = P.tuneSrdkfZGrid;

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
  makeFilter = @(c) SRDKF(plant, Ts, T, netGraph, P.dkfAlpha, P.dkfBeta, P.dkfDelta, ...
                          P.klTolerance, c.z * eye(zDim), 'closed');
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
    'totalRuns', totalTuneRuns, 'filterName', 'SRDKF-Closed', ...
    'zGrid', zGrid, ...
    'bases', struct('z', NaN));
  savedPath = saveRun(sprintf('tuneSRDKFClosed_%s', scenario), P, extras, netGraph, results, struct());

  %% Plotting
  disp("Plotting results")
  plotTuneRun(savedPath);
end
