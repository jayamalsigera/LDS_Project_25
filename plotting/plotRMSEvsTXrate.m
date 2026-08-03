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

  % Per-estimator plotting spec.
  %   Colors grouped by estimator family; line styles separate variants within a family.
  %   Centralized  → blue family    (CKF solid, CRKF dashed)
  %   Consensus    → red/orange     (DSEA-CP solid)
  %   DKF family   → amber          (DKF solid, SDKF-Closed dashed, SDKF-Open dotted)
  %   RDKF family  → purple/magenta (RDKF solid, SRDKF-Closed dashed, SRDKF-Open dotted)
  c_ckf    = [0.00 0.45 0.74];   % MATLAB default blue
  c_crkf   = [0.30 0.65 0.90];   % lighter blue  (distinct from CKF)
  c_dseacp = [0.85 0.33 0.10];   % MATLAB default orange-red
  c_dkf    = [0.93 0.69 0.13];   % MATLAB default amber
  c_rdkf   = [0.49 0.18 0.56];   % MATLAB default purple
  c_rdkfloc = [0.29 0.00 0.51];  % indigo (RDKFLOC — local-tolerance RDKF)
  c_sdkf   = [0.10 0.75 0.65];   % teal  (SDKF family)
  c_srdkf  = [0.47 0.67 0.19];   % green (SRDKF family)
  c_srdkfloc = [0.15 0.40 0.10]; % dark green (SRDKFLOC — local-tolerance SRDKF)

  lw = 1.3;   % base line width

  % spec(label) -> struct with the field prefix, color, line style and width.
  % RMSE data is results.[prefix 'Rmse']; TX rate is results.[prefix 'TxRate'].
  spec = containers.Map('KeyType', 'char', 'ValueType', 'any');
  spec('CKF')          = mkspec('ckf',         c_ckf,    '-',  lw);
  spec('CRKF')         = mkspec('crkf',        c_crkf,   '--', lw);
  spec('DSEA-CP')      = mkspec('dseacp',      c_dseacp, '-',  lw);
  spec('DKF')          = mkspec('dkf',         c_dkf,    '-',  lw);
  spec('RDKF')         = mkspec('rdkf',        c_rdkf,   '-',  lw);
  spec('SDKF-Open')    = mkspec('sdkfOpen',    c_sdkf,   ':',  lw + 0.4);
  spec('SDKF-Closed')  = mkspec('sdkfClosed',  c_sdkf,   '--', lw);
  spec('SRDKF-Open')   = mkspec('srdkfOpen',   c_srdkf,  ':',  lw + 0.4);
  spec('SRDKF-Closed') = mkspec('srdkfClosed', c_srdkf,  '--', lw);
  spec('RDKFLOC')       = mkspec('rdkfloc',       c_rdkfloc,  '-',  lw);
  spec('SRDKFLOC-Open')   = mkspec('srdkflocOpen',   c_srdkfloc, ':',  lw + 0.4);
  spec('SRDKFLOC-Closed') = mkspec('srdkflocClosed', c_srdkfloc, '--', lw);

  % Comparison groups: {name, members}.
  groups = {
    'CKF vs CRKF',                    {'CKF', 'CRKF'};
    'CKF vs DSEA-CP vs DKF',          {'CKF', 'DSEA-CP', 'DKF'};
    'DKF vs RDKF',                    {'DKF', 'RDKF'};
    'DKF vs SDKF-Open vs SRDKF-Open', {'DKF', 'SDKF-Open', 'SRDKF-Open'};
    'DKF vs SDKF-Closed vs SRDKF-Closed', {'DKF', 'SDKF-Closed', 'SRDKF-Closed'};
  };

  % Local-tolerance (Algorithm 2) comparisons: uniform-b robust vs per-node b^i,
  % against the b=0 (DKF) and centralized (CRKF) references. Appended only when
  % the run actually contains the LOC filters, so nominal/older runs are
  % unaffected.
  if isfield(results, 'rdkflocRmse')
    groups = [groups; {
      'DKF vs RDKF vs RDKFLOC',                  {'DKF', 'RDKF', 'RDKFLOC'};
      'DKF vs SRDKF-Open vs SRDKFLOC-Open',      {'DKF', 'SRDKF-Open', 'SRDKFLOC-Open'};
      'DKF vs SRDKF-Closed vs SRDKFLOC-Closed',  {'DKF', 'SRDKF-Closed', 'SRDKFLOC-Closed'};
      'CRKF vs RDKFLOC vs SRDKFLOC-Open',        {'CRKF', 'RDKFLOC', 'SRDKFLOC-Open'};
    }];
  end

  t = (0:T) * Ts;

  for g = 1:size(groups, 1)
    name    = groups{g, 1};
    members = groups{g, 2};

    rmseMembers = presentMembers(results, spec, members, 'Rmse');
    % A single-line TX-rate plot (e.g. DKF alone) is not worth showing.
    txMembers   = presentMembers(results, spec, members, 'TxRate');
    if numel(txMembers) < 2; txMembers = {}; end

    if ~isempty(rmseMembers) && ~isempty(txMembers)
      % Both metrics available: side-by-side subplots.
      figure
      subplot(1, 2, 1);
      drawMetric(t, results, spec, rmseMembers, 'Rmse', name, isLfm);
      subplot(1, 2, 2);
      drawMetric(t, results, spec, txMembers, 'TxRate', name, isLfm);
    else
      if ~isempty(rmseMembers)
        figure; drawMetric(t, results, spec, rmseMembers, 'Rmse', name, isLfm);
      end
      if ~isempty(txMembers)
        figure; drawMetric(t, results, spec, txMembers, 'TxRate', name, isLfm);
      end
    end
  end
end

function s = mkspec(prefix, color, style, lw)
  s = struct('prefix', prefix, 'color', color, 'style', style, 'lw', lw);
end

function present = presentMembers(results, spec, members, suffix)
% Members of a group for which the requested metric field exists.
  present = {};
  for k = 1:numel(members)
    s = spec(members{k});
    if isfield(results, [s.prefix suffix])
      present{end+1} = members{k}; %#ok<AGROW>
    end
  end
end

function drawMetric(t, results, spec, members, suffix, groupName, isLfm)
% Draw one metric ('Rmse' or 'TxRate') into the current axes.
  isRmse = strcmp(suffix, 'Rmse');

  hold on
  for k = 1:numel(members)
    s = spec(members{k});
    y = mean(results.([s.prefix suffix]), 1);
    if isRmse
      semilogy(t, y, s.style, 'Color', s.color, 'LineWidth', s.lw, 'DisplayName', members{k});
    else
      plot(t, y, s.style, 'Color', s.color, 'LineWidth', s.lw, 'DisplayName', members{k});
    end
  end
  hold off

  if isRmse
    metricName = 'RMSE'; ylab = 'RMSE';
    set(gca, 'YScale', 'log');
  else
    metricName = 'TX Rate'; ylab = 'TX Rate';
  end
  title(sprintf('%s vs Time — %s%s', metricName, groupName, optstr(isLfm, ' (LFM data)', '')));
  xlabel('Time (s)'); ylabel(ylab);
  legend('Location', 'northeast'); grid();
end

function s = optstr(cond, a, b)
  if cond; s = a; else; s = b; end
end
