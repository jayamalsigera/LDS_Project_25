function plotTuneRun(varargin)
%PLOTTUNERUN  Re-create figures for one or more saved tune* runs.
%
%   plotTuneRun(path)
%       For a single run: per-sweep RMSE/TX-vs-time figures plus the
%       RMSE-vs-TX-rate tradeoff curve.
%
%   plotTuneRun(path1, path2, ...)
%       Overlay the tradeoff curves of several runs on one axis (e.g. a
%       matched-rate SDKF-Closed vs DKF comparison). The per-sweep time
%       figures are skipped to keep the overlay readable.
%
% Parameter-agnostic: the set of swept parameters is discovered from the
% stored configs; each config carries a 'sweep' field naming its varied
% parameter. When a run stores steady-state RMSE stats (ssRmseMean /
% ssRmseStd) the tradeoff uses them with +/- 1 std error bars; otherwise it
% falls back to the time-averaged mean RMSE with no bars (older runs).

  paths = varargin;
  if isempty(paths)
    error('plotTuneRun:noInput', 'Provide at least one saved tune* .mat path.');
  end

  % A single delta x beta grid run gets its own scatter + heatmap figures
  % (the OAT 'sweep'=field-name convention below does not apply to grids).
  if numel(paths) == 1
    runData = loadRun(paths{1});
    if isfield(runData.extras, 'gridAxisRow') && isfield(runData.extras, 'gridAxisCol')
      plotGridRun(runData);
      return;
    end
    % Per-sweep time-series figures only make sense for a single run.
    plotTimeSeriesFigures(runData);
  end

  % Combined RMSE-vs-TX tradeoff (works for 1..N runs).
  styles = {'-', '--', ':', '-.'};
  figure; hold on;
  allSs = true;
  for k = 1:numel(paths)
    runData = loadRun(paths{k});
    r       = runData.results;
    fn      = runData.extras.filterName;
    configs = r.configs;
    hasSs   = isfield(r, 'ssRmseMean');
    allSs   = allSs && hasSs;
    baseCol = familyColor(fn, k);

    sweeps = unique(cellfun(@(c) string(c.sweep), configs, 'UniformOutput', true));
    multi  = (numel(paths) > 1) || (numel(sweeps) > 1);
    for si = 1:numel(sweeps)
      sw  = char(sweeps(si));
      idx = find(cellfun(@(c) strcmp(c.sweep, sw), configs));

      x = r.meanTxRate(idx);
      if hasSs
        y = r.ssRmseMean(idx); yerr = r.ssRmseStd(idx);
      else
        y = r.meanRmse(idx);   yerr = zeros(numel(idx), 1);
      end
      knob = arrayfun(@(j) configs{j}.(sw), idx);

      % Sort by TX rate so the connecting line runs left-to-right.
      [x, ord] = sort(x(:));
      y = y(ord); yerr = yerr(ord); knob = knob(ord);

      style = styles{mod(si - 1, numel(styles)) + 1};
      if multi, label = sprintf('%s / %s', fn, sw); else, label = fn; end
      errorbar(x, y, yerr, ['o' style], 'Color', baseCol, ...
               'MarkerFaceColor', baseCol, 'LineWidth', 1.3, 'CapSize', 4, ...
               'DisplayName', label);
      for i = 1:numel(x)
        text(x(i), y(i), sprintf('  %s=%.2g', sw, knob(i)), ...
             'Color', baseCol, 'FontSize', 7, 'VerticalAlignment', 'bottom');
      end
    end
  end
  hold off;
  xlabel('Mean Transmission Rate');
  if allSs, ylabel('Steady-state RMSE'); else, ylabel('Mean RMSE'); end
  title('RMSE vs Transmission Rate Tradeoff');
  legend('Location', 'best'); grid on;
end

function plotTimeSeriesFigures(runData)
% Per-sweep RMSE and TX rate vs time (one figure per swept parameter).
  params  = runData.params;
  extras  = runData.extras;
  results = runData.results;

  T  = params.T;
  Ts = params.Ts;
  t  = (0:T) * Ts;

  configs    = results.configs;
  rmseCurves = results.rmseCurves;
  txCurves   = results.txCurves;
  nConfigs   = numel(configs);

  filterName = extras.filterName;
  bases      = extras.bases;

  sweepNames = unique(arrayfun(@(k) string(configs{k}.sweep), 1:nConfigs));
  for k = 1:numel(sweepNames)
    sweepName  = char(sweepNames(k));
    fixedNames = setdiff(fieldnames(bases), {sweepName});
    fixedPairs = cell(1, 2 * numel(fixedNames));
    for j = 1:numel(fixedNames)
      fixedPairs{2*j-1} = fixedNames{j};
      fixedPairs{2*j}   = bases.(fixedNames{j});
    end
    plotSweep(t, configs, rmseCurves, txCurves, sweepName, fixedPairs, filterName);
  end
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

function plotGridRun(runData)
% Two figures for a full-factorial (row x col) grid run:
%   1. RMSE-vs-TX scatter with the Pareto frontier highlighted;
%   2. RMSE and TX heatmaps over the (row, col) grid.
  r      = runData.results;
  extras = runData.extras;
  fn     = extras.filterName;
  col    = familyColor(fn, 1);

  rowAx = extras.gridAxisRow;   % inner loop -> reshape rows
  colAx = extras.gridAxisCol;   % outer loop -> reshape cols
  nRow  = numel(rowAx.values);
  nCol  = numel(colAx.values);

  x    = r.meanTxRate(:);
  hasSs = isfield(r, 'ssRmseMean');
  if hasSs
    y = r.ssRmseMean(:); yerr = r.ssRmseStd(:); yLab = 'Steady-state RMSE';
  else
    y = r.meanRmse(:);   yerr = zeros(size(y)); yLab = 'Mean RMSE';
  end
  if isfield(r, 'onFront'), front = logical(r.onFront(:)); else, front = paretoFront(x, y); end

  configs = r.configs;
  rowVal = arrayfun(@(k) configs{k}.(rowAx.name), 1:numel(configs))';
  colVal = arrayfun(@(k) configs{k}.(colAx.name), 1:numel(configs))';

  %% Figure 1: tradeoff scatter + Pareto frontier
  figure; hold on;
  errorbar(x, y, yerr, 'o', 'Color', col, 'MarkerFaceColor', col, ...
           'LineStyle', 'none', 'CapSize', 4, 'DisplayName', [fn ' configs']);
  [xf, ord] = sort(x(front));
  yf = y(front); yf = yf(ord);
  plot(xf, yf, '-', 'Color', col, 'LineWidth', 1.6, 'DisplayName', 'Pareto frontier');
  for i = 1:numel(x)
    text(x(i), y(i), sprintf('  %s=%.2g, %s=%.2g', colAx.name, colVal(i), rowAx.name, rowVal(i)), ...
         'Color', col, 'FontSize', 7, 'VerticalAlignment', 'bottom');
  end
  hold off;
  xlabel('Mean Transmission Rate'); ylabel(yLab);
  title(sprintf('%s: RMSE vs TX tradeoff (%s x %s grid)', fn, rowAx.name, colAx.name));
  legend('Location', 'best'); grid on;

  %% Figure 2: RMSE and TX heatmaps over the grid
  % Config order is col outer, row inner -> column-major reshape(v, nRow, nCol).
  Mrmse = reshape(y, nRow, nCol);
  Mtx   = reshape(x, nRow, nCol);
  figure;
  gridHeatmap(subplot(1, 2, 1), Mrmse, colAx, rowAx, sprintf('%s %s', fn, yLab));
  gridHeatmap(subplot(1, 2, 2), Mtx,   colAx, rowAx, sprintf('%s Mean TX Rate', fn));
end

function gridHeatmap(ax, M, colAx, rowAx, ttl)
  imagesc(ax, M);
  set(ax, 'XTick', 1:numel(colAx.values), 'XTickLabel', colAx.values, ...
          'YTick', 1:numel(rowAx.values), 'YTickLabel', rowAx.values, ...
          'YDir', 'normal');
  xlabel(ax, colAx.name); ylabel(ax, rowAx.name);
  title(ax, ttl); colorbar(ax);
  for ci = 1:numel(colAx.values)
    for ri = 1:numel(rowAx.values)
      text(ax, ci, ri, sprintf('%.2g', M(ri, ci)), ...
           'HorizontalAlignment', 'center', 'FontSize', 7, 'Color', [0 0 0]);
    end
  end
end

function col = familyColor(filterName, fallbackIdx)
% Repo family colours (see plotSST2dRun); fallback palette otherwise.
  switch filterName
    case 'DKF'
      col = [0.93 0.69 0.13];   % amber
    case {'SDKF-Closed', 'SDKF-Open'}
      col = [0.10 0.75 0.65];   % teal
    case {'SRDKF-Closed', 'SRDKF-Open'}
      col = [0.47 0.67 0.19];   % green
    case 'RDKF'
      col = [0.49 0.18 0.56];   % purple
    otherwise
      pal = lines(max(fallbackIdx, 7));
      col = pal(fallbackIdx, :);
  end
end
