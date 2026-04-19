function plotTuneRun(runData)
%PLOTTUNERUN  Re-create figures for a saved tuneDKF/tuneRDKF run.
%
% Handles both DKF (alpha/beta/delta) and RDKF (alpha/beta/delta/b) sweeps;
% the extra b-parameter panel is emitted when the stored configs carry a
% 'b' field.

  params   = runData.params;
  results  = runData.results;
  extras   = runData.extras;

  T  = params.T;
  Ts = params.Ts;
  t  = (0:T) * Ts;

  configs    = results.configs;
  rmseCurves = results.rmseCurves;
  txCurves   = results.txCurves;
  meanRmse   = results.meanRmse;
  meanTxRate = results.meanTxRate;
  nConfigs   = numel(configs);

  hasB = isfield(configs{1}, 'b');
  if hasB
    filterName = 'RDKF';
  else
    filterName = 'DKF';
  end

  if hasB
    plotSweep(t, configs, rmseCurves, txCurves, 'alpha', ...
              {'beta', extras.betaBase, 'delta', extras.deltaBase, 'b', extras.bBase}, filterName);
    plotSweep(t, configs, rmseCurves, txCurves, 'beta', ...
              {'alpha', extras.alphaBase, 'delta', extras.deltaBase, 'b', extras.bBase}, filterName);
    plotSweep(t, configs, rmseCurves, txCurves, 'delta', ...
              {'alpha', extras.alphaBase, 'beta', extras.betaBase, 'b', extras.bBase}, filterName);
    plotSweep(t, configs, rmseCurves, txCurves, 'b', ...
              {'alpha', extras.alphaBase, 'beta', extras.betaBase, 'delta', extras.deltaBase}, filterName);
    sweepColors = struct('alpha', [0.85 0.33 0.10], ...
                         'beta',  [0.00 0.45 0.74], ...
                         'delta', [0.47 0.67 0.19], ...
                         'b',     [0.49 0.18 0.56]);
  else
    plotSweep(t, configs, rmseCurves, txCurves, 'alpha', ...
              {'beta', extras.betaBase, 'delta', extras.deltaBase}, filterName);
    plotSweep(t, configs, rmseCurves, txCurves, 'beta', ...
              {'alpha', extras.alphaBase, 'delta', extras.deltaBase}, filterName);
    plotSweep(t, configs, rmseCurves, txCurves, 'delta', ...
              {'alpha', extras.alphaBase, 'beta', extras.betaBase}, filterName);
    sweepColors = struct('alpha', [0.85 0.33 0.10], ...
                         'beta',  [0.00 0.45 0.74], ...
                         'delta', [0.47 0.67 0.19]);
  end

  % RMSE vs TX rate tradeoff scatter
  figure;
  hold on;
  for i = 1:nConfigs
    c = configs{i};
    col = sweepColors.(c.sweep);
    scatter(meanTxRate(i), meanRmse(i), 60, col, 'filled', ...
            'DisplayName', sprintf('%s=%.2g', c.sweep, c.(c.sweep)));
    text(meanTxRate(i), meanRmse(i), ...
         sprintf(' %s=%.2g', c.sweep, c.(c.sweep)), 'FontSize', 8);
  end
  hold off;
  xlabel('Mean Transmission Rate');
  ylabel('Mean RMSE');
  title(sprintf('%s: RMSE vs Transmission Rate Tradeoff', filterName));
  grid on;
end

function plotSweep(t, configs, rmseCurves, txCurves, param, fixedPairs, filterName)
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
  title(sprintf('%s RMSE sweep over %s (others fixed: %s)', filterName, param, fixedStr));
  legend(); grid on;

  subplot(2, 1, 2);
  hold on;
  for i = idx
    v = configs{i}.(param);
    plot(t, txCurves(i, :), 'DisplayName', sprintf('%s=%.2g', param, v));
  end
  hold off;
  xlabel('Time (s)'); ylabel('TX Rate');
  title(sprintf('%s Transmission Rate sweep over %s', filterName, param));
  legend(); grid on;
end
