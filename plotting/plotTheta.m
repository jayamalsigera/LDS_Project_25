function plotTheta(path)
% plotTheta Re-create the Risk Sensitivity Parameter (theta) vs Time figure
% for a saved estimateSST run.
%
% Usage:
%   plotTheta('results/some_run.mat')
%   plotTheta() % Loads the most recent run

  if nargin < 1
    path = '';
  end
  runData = loadRun(path);

  params   = runData.params;
  results  = runData.results;

  T  = params.T;
  Ts = params.Ts;
  isLfm = contains(runData.script, 'Lfm');

  % Spec for robust estimators only
  spec = {
    'CRKF',             'crkf',             [0.30 0.65 0.90], '--'
    'RDKF',             'rdkf',             [0.49 0.18 0.56], '-'
    'RDKFLOC',          'rdkfloc',          [0.29 0.00 0.51], '-'
    'SRDKF-Open',       'srdkfOpen',        [0.47 0.67 0.19], ':'
    'SRDKF-Closed',     'srdkfClosed',      [0.47 0.67 0.19], '--'
    'SRDKFLOC-Open',    'srdkflocOpen',     [0.15 0.40 0.10], ':'
    'SRDKFLOC-Closed',  'srdkflocClosed',   [0.15 0.40 0.10], '--'
  };

  t = (0:T) * Ts;

  figure;
  hold on;

  for k = 1:size(spec, 1)
      label  = spec{k, 1};
      prefix = spec{k, 2};
      color  = spec{k, 3};
      style  = spec{k, 4};

      field = [prefix 'Theta'];
      if ~isfield(results, field)
          continue;
      end

      data = results.(field);

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

      plot(t, y, style, 'Color', color, 'LineWidth', 1.5, 'DisplayName', label);
  end

  hold off;
  xlabel('Time (s)');
  ylabel('Risk sensitivity parameter \theta');
  if isLfm
      title('Risk Sensitivity \theta vs Time (LFM data)');
  else
      title('Risk Sensitivity \theta vs Time');
  end
  legend('Location', 'northeast');
  grid on;
end
