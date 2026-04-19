function plotTuneRun(path)
%PLOTTUNERUN  Re-create figures for a saved tune* run.
%
% Parameter-agnostic: the set of swept parameters is discovered from the
% stored configs. Each config must carry a 'sweep' field naming the
% parameter it varies; extras must carry 'filterName' (string) and
% 'bases' (struct of param -> baseline value).

  runData = loadRun(path);

  params   = runData.params;
  extras   = runData.extras;
  results  = runData.results;

  T  = params.T;
  Ts = params.Ts;
  t  = (0:T) * Ts;

  configs    = results.configs;
  rmseCurves = results.rmseCurves;
  txCurves   = results.txCurves;
  meanRmse   = results.meanRmse;
  meanTxRate = results.meanTxRate;
  nConfigs   = numel(configs);

  filterName = extras.filterName;
  bases      = extras.bases;

  sweepNames = unique(arrayfun(@(k) string(configs{k}.sweep), 1:nConfigs));
  palette = lines(numel(sweepNames));
  sweepColors = struct();
  for k = 1:numel(sweepNames)
    sweepColors.(char(sweepNames(k))) = palette(k, :);
  end

  for k = 1:numel(sweepNames)
    sweepName = char(sweepNames(k));
    fixedNames = setdiff(fieldnames(bases), {sweepName});
    fixedPairs = cell(1, 2 * numel(fixedNames));
    for j = 1:numel(fixedNames)
      fixedPairs{2*j-1} = fixedNames{j};
      fixedPairs{2*j}   = bases.(fixedNames{j});
    end
    plotSweep(t, configs, rmseCurves, txCurves, sweepName, fixedPairs, filterName);
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
