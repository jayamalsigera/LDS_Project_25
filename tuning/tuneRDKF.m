%% RDKF Hyperparameter Tuning (delta x beta grid)
%
%   tuneRDKF()
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
% layer's knob, not a trigger knob, and is left at its baseline here. Grids
% live in sst3dParams.m.
%
function tuneRDKF()

  rng(42);

  %% Fixed parameters
  P = collectParams();
  Ts = P.Ts; T = P.T;
  totalTuneRuns = P.totalTuneRuns;

  %% Hyperparameter grid (delta x beta; alpha and b fixed)
  betaGrid   = P.tuneRdkfBetaGrid;
  deltaGrid  = P.tuneRdkfDeltaGrid;
  alphaFixed = P.tuneRdkfAlphaFixed;
  bFixed     = P.tuneRdkfBFixed;

  %% Network and Plant
  disp("Creating Network")
  netGraph = createSpatialNetwork(P.nodeCount, P.sensorCount, P.maxLength);

  disp("Simulating target dynamics")
  plant = SingleTarget3dModel(P.Ts, P.sensorCount, P.noiseScale, P.T);

  %% Pre-generate trajectories so all configurations see the same data
  disp("Pre-generating Monte Carlo trajectories")
  samples = cell(totalTuneRuns, 1);
  for run = 1:totalTuneRuns
    samples{run} = plant.simulate(P.x0);
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

  disp("Running simulations...")
  makeFilter = @(c) RDKF(plant, Ts, T, netGraph, c.alpha, c.beta, c.delta, c.b);
  [meanRmse, finalRmse, meanTxRate, ssRmseMean, ssRmseStd, rmseCurves, txCurves] = ...
      evalConfigsMC(makeFilter, configs, samples, P.x0_hat, P.P0, T);

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
  savedPath = saveRun(mfilename, P, extras, netGraph, results);

  %% Plotting
  disp("Plotting results")
  plotTuneRun(savedPath);
end
