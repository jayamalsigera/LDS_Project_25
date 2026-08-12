function plotRMSEvsTXrate(path)
%PLOTRMSEVSTXRATE  Re-create RMSE / TX-rate figures for a saved estimateSST run.
%
% Works for the nominal-plant script (estimateSST3d) and its least-favorable
% model counterpart (estimateSST3dLfm); the LFM branch is taken when
% runData.samples has a nomSample field.

  runData = loadRun(path);

  params   = runData.params;
  netGraph = runData.netGraph;
  results  = runData.results;
  samples  = runData.samples;

  T  = params.T;
  Ts = params.Ts;
  isLfm = isfield(samples, 'nomSample');

  % Network layout
  plotNetwork(netGraph, params.maxLength);

  % Per-estimator plotting spec: {label, results field prefix, color, line style, hasTxRate}.
  %   Colors grouped by estimator family; line styles separate variants within a family.
  %   Centralized  → blue family    (CKF solid, CRKF dashed)
  %   Consensus    → red/orange     (DSEA-CP solid)
  %   DKF family   → amber          (DKF solid, SDKF-Closed dashed, SDKF-Open dotted)
  %   RDKF family  → purple/magenta (RDKF solid, SRDKF-Closed dashed, SRDKF-Open dotted)
  %   LOC variants → indigo / dark green (per-node local tolerances b^i)
  % RMSE data is results.[prefix 'Rmse']; TX rate is results.[prefix 'TxRate'],
  % which only the event-triggered / DKF-family filters record.
  spec = {
    'CKF',              'ckf',              [0.00 0.45 0.74], '-',  false
    'CRKF',             'crkf',             [0.30 0.65 0.90], '--', false
    'DSEA-CP',          'dseacp',           [0.85 0.33 0.10], '-',  false
    'DKF',              'dkf',              [0.93 0.69 0.13], '-',  true
    'RDKF',             'rdkf',             [0.49 0.18 0.56], '-',  true
    'RDKFLOC',          'rdkfloc',          [0.29 0.00 0.51], '-',  true
    'SDKF-Open',        'sdkfOpen',         [0.10 0.75 0.65], ':',  true
    'SDKF-Closed',      'sdkfClosed',       [0.10 0.75 0.65], '--', true
    'SRDKF-Open',       'srdkfOpen',        [0.47 0.67 0.19], ':',  true
    'SRDKF-Closed',     'srdkfClosed',      [0.47 0.67 0.19], '--', true
    'SRDKFLOC-Open',    'srdkflocOpen',     [0.15 0.40 0.10], ':',  true
    'SRDKFLOC-Closed',  'srdkflocClosed',   [0.15 0.40 0.10], '--', true
  };

  % Comparison groups: {name, members}.
  groups = {
    'CKF vs CRKF',                             {'CKF', 'CRKF'};
    'CKF vs DSEA-CP vs DKF',                   {'CKF', 'DSEA-CP', 'DKF'};
    'DKF vs SDKF-Open vs SRDKF-Open',          {'DKF', 'SDKF-Open', 'SRDKF-Open'};
    'DKF vs SDKF-Closed vs SRDKF-Closed',      {'DKF', 'SDKF-Closed', 'SRDKF-Closed'};
    'DKF vs RDKF vs RDKFLOC',                  {'DKF', 'RDKF', 'RDKFLOC'};
    'DKF vs SRDKF-Open vs SRDKFLOC-Open',      {'DKF', 'SRDKF-Open', 'SRDKFLOC-Open'};
    'DKF vs SRDKF-Closed vs SRDKFLOC-Closed',  {'DKF', 'SRDKF-Closed', 'SRDKFLOC-Closed'};
    'CRKF vs RDKFLOC vs SRDKFLOC-Open',        {'CRKF', 'RDKFLOC', 'SRDKFLOC-Open'};
  };

  t = (0:T) * Ts;

  for g = 1:size(groups, 1)
    name    = groups{g, 1};
    members = groups{g, 2};

    % A single-line TX-rate plot (e.g. DKF alone) is not worth showing.
    txMembers = members(cellfun(@(m) spec{strcmp(spec(:, 1), m), 5}, members));
    if numel(txMembers) < 2
      figure; drawMetric(t, results, spec, members, 'Rmse', name, isLfm);
    else
      figure
      subplot(1, 2, 1);
      drawMetric(t, results, spec, members, 'Rmse', name, isLfm);
      subplot(1, 2, 2);
      drawMetric(t, results, spec, txMembers, 'TxRate', name, isLfm);
    end
  end
end

function drawMetric(t, results, spec, members, suffix, groupName, isLfm)
% Draw one metric ('Rmse' or 'TxRate') into the current axes.
  isRmse = strcmp(suffix, 'Rmse');

  hold on
  for k = 1:numel(members)
    row    = strcmp(spec(:, 1), members{k});
    prefix = spec{row, 2};
    color  = spec{row, 3};
    style  = spec{row, 4};

    y = mean(results.([prefix suffix]), 1);
    if isRmse
      semilogy(t, y, style, 'Color', color, 'LineWidth', 1.3, 'DisplayName', members{k});
    else
      plot(t, y, style, 'Color', color, 'LineWidth', 1.3, 'DisplayName', members{k});
    end
  end
  hold off

  if isRmse
    metricName = 'RMSE';
    set(gca, 'YScale', 'log');
  else
    metricName = 'TX Rate';
  end
  if isLfm
    titleSuffix = ' (LFM data)';
  else
    titleSuffix = '';
  end
  title(sprintf('%s vs Time — %s%s', metricName, groupName, titleSuffix));
  xlabel('Time (s)'); ylabel(metricName);
  legend('Location', 'northeast'); grid();
end
