function plotTheta(path)
% plotTheta Re-create the Risk Sensitivity Parameter (theta vs theta_bar) figures
% for each robust estimator from a saved estimateSST run.
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

  % Spec for the 6 robust distributed estimators: {Label, results field prefix, color}
  spec = {
    'RDKF',             'rdkf',             [0.49 0.18 0.56]
    'RDKFLOC',          'rdkfloc',          [0.29 0.00 0.51]
    'SRDKF',            'srdkf',            [0.47 0.67 0.19]
    'SRDKFLOC',         'srdkfloc',         [0.15 0.40 0.10]
  };

  t = (0:T) * Ts;

  if ~exist('results/figures', 'dir')
    mkdir('results/figures');
  end

  % Extract average curves for theta and theta_bar to determine unified y-limits
  numEst = size(spec, 1);
  yTheta = cell(numEst, 1);
  yThetaBar = cell(numEst, 1);
  globalMin = Inf;
  globalMax = -Inf;

  for k = 1:numEst
    prefix = spec{k, 2};

    fieldTheta = [prefix 'Theta'];
    if isfield(results, fieldTheta)
      y = computeAverage(results.(fieldTheta));
      yTheta{k} = y;
      posVals = y(y > 0 & isfinite(y));
      if ~isempty(posVals)
        globalMin = min(globalMin, min(posVals));
        globalMax = max(globalMax, max(posVals));
      end
    end

    fieldThetaBar = [prefix 'ThetaBar'];
    if isfield(results, fieldThetaBar)
      yBar = computeAverage(results.(fieldThetaBar));
      yThetaBar{k} = yBar;
      posValsBar = yBar(yBar > 0 & isfinite(yBar));
      if ~isempty(posValsBar)
        globalMin = min(globalMin, min(posValsBar));
        globalMax = max(globalMax, max(posValsBar));
      end
    end
  end

  % Define shared y-limits if valid data is present
  if isfinite(globalMin) && isfinite(globalMax) && globalMin < globalMax
    yLimits = [globalMin * 0.9, globalMax * 1.1];
  else
    yLimits = [];
  end

  % Generate individual comparison plots for each estimator
  for k = 1:numEst
    label  = spec{k, 1};
    prefix = spec{k, 2};
    color  = spec{k, 3};

    yT  = yTheta{k};
    yTB = yThetaBar{k};

    if isempty(yT) && isempty(yTB)
      continue;
    end

    figure;
    hold on;

    if ~isempty(yT)
      semilogy(t(2:end), yT(2:end), '-', 'Color', color, 'LineWidth', 1.5, ...
        'DisplayName', '$\theta$ (Local Prior)');
    end
    if ~isempty(yTB)
      semilogy(t(2:end), yTB(2:end), '--', 'Color', color * 0.7, 'LineWidth', 1.5, ...
        'DisplayName', '$\bar{\theta}$ (Global Prior)');
    end

    hold off;
    set(gca, 'YScale', 'log');
    xlabel('Time (s)');
    ylabel('Risk sensitivity parameter');
    labelEsc = strrep(label, '_', '\_');
    if isLfm
      title(sprintf('Risk Sensitivity $\\theta$ vs $\\bar{\\theta}$ - %s (LFM data)', labelEsc), 'Interpreter', 'latex');
    else
      title(sprintf('Risk Sensitivity $\\theta$ vs $\\bar{\\theta}$ - %s', labelEsc), 'Interpreter', 'latex');
    end

    if ~isempty(yLimits)
      ylim(yLimits);
    end

    legend('Location', 'northeast', 'Interpreter', 'latex');
    grid on;

    cleanName = strrep(label, ' ', '_');
    cleanName = strrep(cleanName, '-', '_');
    if isLfm
      exportgraphics(gcf, sprintf('results/figures/theta_comparison_%s_lfm.pdf', cleanName), 'ContentType', 'vector');
    else
      exportgraphics(gcf, sprintf('results/figures/theta_comparison_%s.pdf', cleanName), 'ContentType', 'vector');
    end
  end
end

function y = computeAverage(data)
  % Compute run- and node-averaged curve
  if ndims(data) == 2 && size(data, 1) == 1
    y = data;
  elseif ndims(data) == 2
    y = mean(data, 1);
  else
    % Distributed data: (totalRuns, nodeCount, T+1)
    avgRuns = squeeze(mean(data, 1));
    if isvector(avgRuns)
      y = avgRuns;
    else
      y = mean(avgRuns, 1);
    end
  end
end
