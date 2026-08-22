function plotTheta(path)
% plotTheta Re-create the Risk Sensitivity Parameter (theta and theta_bar) vs Time
% figures for a saved estimateSST run.
%
% Usage:
%   plotTheta('results/some_run.mat')
%   plotTheta() % Loads the most recent run

  if nargin < 1
    path = '';
  end
  runData = loadRun(path, 'estimateSST*.mat');

  params   = runData.params;
  results  = runData.results;

  T  = params.T;
  Ts = params.Ts;
  isLfm = contains(runData.script, 'Lfm');

  % Spec for robust estimators only
  spec = {
    % 'CRKF',             'crkf',             [0.30 0.65 0.90], '--'
    'RDKF',             'rdkf',             [0.49 0.18 0.56], '-'
    'RDKFLOC',          'rdkfloc',          [0.29 0.00 0.51], '-'
    'SRDKF-Open',       'srdkfOpen',        [0.47 0.67 0.19], ':'
    'SRDKF-Closed',     'srdkfClosed',      [0.47 0.67 0.19], '--'
    'SRDKFLOC-Open',    'srdkflocOpen',     [0.15 0.40 0.10], ':'
    'SRDKFLOC-Closed',  'srdkflocClosed',   [0.15 0.40 0.10], '--'
  };

  t = (0:T) * Ts;

  if ~exist('results/figures', 'dir')
    mkdir('results/figures');
  end

  % Plot theta (local priors)
  plotMetric(t, results, spec, 'Theta', '\theta', 'theta_vs_time', isLfm);

  % Plot theta_bar (global/no-transmit priors)
  plotMetric(t, results, spec, 'ThetaBar', '\bar{\theta}', 'theta_bar_vs_time', isLfm);
end

function plotMetric(t, results, spec, suffix, symbolStr, fileBaseName, isLfm)
  figure;
  hold on;
  hasData = false;

  for k = 1:size(spec, 1)
      label  = spec{k, 1};
      prefix = spec{k, 2};
      color  = spec{k, 3};
      style  = spec{k, 4};

      field = [prefix suffix];
      if ~isfield(results, field)
          continue;
      end

      data = results.(field);
      hasData = true;

      if strcmp(prefix, 'crkf')
          % CRKF data is (totalRuns, T+1)
          y = mean(data, 1);
      else
          % Distributed data is (totalRuns, nodeCount, T+1)
          % 1. Average across runs -> size: (1, nodeCount, T+1) -> squeeze to (nodeCount, T+1)
          avgRuns = squeeze(mean(data, 1));

          % 2. Average across nodes
          % If nodeCount is 1, avgRuns is 1D (T+1), so check size
          if isvector(avgRuns)
              y = avgRuns;
          else
              y = mean(avgRuns, 1);
          end
      end

      semilogy(t, y, style, 'Color', color, 'LineWidth', 1.5, 'DisplayName', label);
  end

  hold off;
  if ~hasData
      close(gcf);
      return;
  end

  xlabel('Time (s)');
  ylabel(sprintf('Risk sensitivity parameter %s', symbolStr));
  if isLfm
      title(sprintf('Risk Sensitivity %s vs Time (LFM data)', symbolStr));
  else
      title(sprintf('Risk Sensitivity %s vs Time', symbolStr));
  end
  legend('Location', 'northeast');
  grid on;

  if isLfm
    exportgraphics(gcf, sprintf('results/figures/%s_lfm.pdf', fileBaseName), 'ContentType', 'vector');
  else
    exportgraphics(gcf, sprintf('results/figures/%s.pdf', fileBaseName), 'ContentType', 'vector');
  end
end
