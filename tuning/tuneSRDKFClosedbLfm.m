function tuneSRDKFClosedbLfm()
%% SRDKF (closed-loop trigger) robust-tolerance (b) tuning on LEAST-FAVORABLE data
%
%   tuneSRDKFClosedbLfm()
%
% Sweeps the KL-divergence tolerance b (a.k.a. klTolerance) -- the robust
% layer's knob that sizes the per-node deflation Psi = Omega_pred - theta*I --
% with the stochastic-trigger weight (Z = errorNormWeightsClosed) and the
% other trigger knobs held at the sst3dParams.m baseline. This isolates the
% robustness knob that tuneSRDKFClosed leaves fixed.
%
% Two differences from tuneSRDKFClosed (which sweeps z):
%   * the swept axis is b (grid tuneBGrid), not the error-norm weight scale z;
%   * trajectories are drawn from the LEAST-FAVORABLE model (as in
%     estimateSST3dLfm.m), not the nominal plant -- robustness only earns its
%     keep under model mismatch, so on nominal data the optimal b is ~0.
%
% Headless/server run: NO plotting. Move the saved .mat off the cluster and
% plot locally with plotTuneRun(savedPath).

  rng(42);

  %% Fixed parameters
  P = collectParams();
  Ts = P.Ts; T = P.T;
  totalTuneRuns = P.totalTuneRuns;

  %% Robust-tolerance grid
  % b = 0 is the anchor (theta = 0 => Psi = Omega_pred => robust layer disabled);
  % the rest is a log grid through and just past the shipped 0.05.
  bGrid = P.tuneBGrid;

  %% Network and Plant
  disp("Creating Network")
  netGraph = createSpatialNetwork(P.nodeCount, P.sensorCount, P.maxLength);

  disp("Building plant and least-favorable model")
  plant = SingleTarget3dModel(P.Ts, P.sensorCount, P.noiseScale, P.T);
  lfm   = LeastFavorableModel(plant, P.P0, P.lfmKlTolerance, T);

  %% Pre-generate LFM trajectories so all configurations see the same data
  disp("Pre-generating least-favorable Monte Carlo trajectories")
  samples = cell(totalTuneRuns, 1);
  for run = 1:totalTuneRuns
    % Keep only the (X, Y) trajectories. lfm.simulate returns the whole
    % LeastFavorableModel object, whose per-step sweep arrays (notably L, an
    % m-by-m-by-(T+1) factor ~ hundreds of MB) would otherwise be stored 50x and
    % broadcast to every parfor worker -- blowing past node memory.
    s = lfm.simulate(P.x0);
    samples{run} = struct('X', s.X, 'Y', s.Y);
  end

  %% Sweep
  nConfigs = numel(bGrid);
  configs = cell(nConfigs, 1);
  for k = 1:nConfigs
    configs{k} = struct('b', bGrid(k), 'sweep', 'b');
  end

  disp("Running simulations...")
  makeFilter = @(c) SRDKF(plant, Ts, T, netGraph, P.dkfAlpha, P.dkfBeta, P.dkfDelta, ...
                          c.b, P.errorNormWeightsClosed, 'closed');
  [meanRmse, finalRmse, meanTxRate, ssRmseMean, ssRmseStd, rmseCurves, txCurves] = ...
      evalConfigsMC(makeFilter, configs, samples, P.x0_hat, P.P0, T);

  %% Results table (sorted by b so the response reads top-to-bottom)
  bCol = arrayfun(@(k) configs{k}.b, 1:nConfigs)';
  resultsTable = table(bCol, meanRmse, finalRmse, meanTxRate, ssRmseMean, ssRmseStd, ...
                       'VariableNames', {'b', 'MeanRMSE', 'FinalRMSE', 'MeanTxRate', ...
                                         'SsRMSE', 'SsRMSEStd'});
  resultsTable = sortrows(resultsTable, 'b');
  disp(resultsTable);

  [bestSs, jb] = min(ssRmseMean);
  fprintf('Best b = %.4g : SsRMSE = %.2f, TX = %.3f  (b=0 anchor SsRMSE = %.2f)\n', ...
          bCol(jb), bestSs, meanTxRate(jb), ssRmseMean(bCol == 0));

  %% Save run (no plotting -- headless)
  results = struct( ...
    'rmseCurves', rmseCurves, 'txCurves', txCurves, ...
    'meanRmse', meanRmse, 'finalRmse', finalRmse, 'meanTxRate', meanTxRate, ...
    'ssRmseMean', ssRmseMean, 'ssRmseStd', ssRmseStd, ...
    'configs', {configs}, 'resultsTable', resultsTable);
  extras = struct( ...
    'totalRuns', totalTuneRuns, 'filterName', 'SRDKF-Closed', ...
    'bGrid', bGrid, 'bases', struct('b', NaN));
  savedPath = saveRun(mfilename, P, extras, netGraph, results, struct());
  fprintf('Plot locally with: plotTuneRun(''%s'')\n', savedPath);
end
