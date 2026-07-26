%% RDKF vs DKF on 3D LFM data — the robust-vs-nominal comparison, swept over b
%
% Focused companion to estimateSST3dLfm.m. It carries ONLY the filters needed for
% the robust-vs-nominal question — CKF, CRKF, DKF, RDKF, RDKFLOC — and sweeps the
% KL tolerance b across the feasibility boundary of the event trigger.
%
% Why a separate script (see docs/rdkf-tx-saturation-analysis.md §6.11):
%
%  1. It runs at |S| = 20 sensor nodes + 80 communication nodes, which is the
%     configuration Ghion & Zorzi (2023) Section 6 actually uses. The shared
%     sst3dParams.m pins sensorCount = nodeCount = 100 because the STOCHASTIC
%     trigger implementations (SDKF/SRDKF/SRDKFLOC) cannot handle non-sensor nodes.
%     Those filters are deliberately absent here, so the |S| = 20 override is safe;
%     estimateSST3dLfm.m must keep using |S| = 100 until that is resolved
%     (checklist item 0).
%
%  2. b is SHARED by CRKF and RDKF on purpose. At a given b the centralized robust
%     filter defends exactly the same KL ball as the distributed one, which makes
%     "does the distributed filter capture the same robustness benefit?" a direct
%     comparison rather than an inference.
%
%  3. The b grid straddles the trigger's feasibility bound. Condition (F) of §6.11,
%     lambda_bar(b) <= 1 + beta, gives b_max = 0.0088 at beta = 0.2: below it eq. (9)
%     can fire, above it the trigger is dead and RDKF pays the deflation cost with no
%     bandwidth saving. The grid brackets that value so the crossing is visible.
%
% The Monte Carlo is PAIRED: run r draws its trajectory after rng(trajSeedBase + r),
% so every b and every filter sees bit-identical (X, Y). The margins of interest are
% a few percent against an across-run std of order 0.1, which only a paired
% comparison can resolve.
%
% One result file is written per b, following the usual saveRun naming so each file
% is self-describing (..._N100s20_b0.001_blfm0.05_runs200_...). NOTE: the params
% struct is taken from collectParams and then OVERRIDDEN, because collectParams
% re-sources sst3dParams.m and would otherwise record the unoverridden values --
% i.e. a mislabelled result file.

clear;
close all;
clc;

%% Parameters
sst3dParams;

% --- overrides for this script (see the header) -------------------------------
sensorCount = 20;          % paper's Section 6: 20 sensors + 80 relays
dkfDelta = 0.5;            % paper's value; also DKF's best in the selection sweep
bGrid = [0.001 0.005 0.01 0.05];   % straddles b_max = 0.0088 from condition (F)
lfmKlTolerance = 0.05;     % data mismatch radius, held fixed
trajSeedBase = 1000;       % per-run trajectory seed => pairing across b

fprintf(['Config: N=%d, S=%d, T=%d, runs=%d, alpha=%g, beta=%g, delta=%g\n' ...
         'LFM data radius %g; b grid %s (b_max from (F) = %.4f)\n'], ...
        nodeCount, sensorCount, T, totalRuns, dkfAlpha, dkfBeta, dkfDelta, ...
        lfmKlTolerance, mat2str(bGrid), 0.0088);

%% Network Definition
disp("Creating Network")
rng(42);
netGraph = createSpatialNetwork(nodeCount, sensorCount, maxLength);
assertConnected(netGraph);

%% Plant
disp("Simulating target dynamics")
rng(43);
plant = SingleTarget3dModel(Ts, sensorCount, noiseScale, T);
assertStabilizable(plant.A, plant.B);
assertDetectable(plant.A, plant.C);

senBlock = plant.p / sensorCount;
for i = 1:sensorCount
  idx = senBlock * (i - 1) + (1:senBlock);
  assertLocallyObservable(plant.A, plant.C(idx, :), i);
end

%% Least-favorable data generator (independent of the filters' b)
disp("Precomputing least-favorable model")
lfm = LeastFavorableModel(plant, P0, lfmKlTolerance, T);

%% Monte Carlo, one pass per b
w = round(T / 2) + 1 : T + 1;      % steady-state window
savedPaths = strings(numel(bGrid), 1);
summary = zeros(numel(bGrid), 8);

for ib = 1:numel(bGrid)
  b = bGrid(ib);
  fprintf('\n===== b = %g (%d of %d) =====\n', b, ib, numel(bGrid));

  % Serial construction: RDKFLOC's per-node b^i solve runs once, not per run.
  crkf    = CRKF(plant, Ts, T, b);
  dkf     = DKF(plant, Ts, T, netGraph, dkfAlpha, dkfBeta, dkfDelta);
  rdkf    = RDKF(plant, Ts, T, netGraph, dkfAlpha, dkfBeta, dkfDelta, b);
  rdkfloc = RDKFLOC(plant, Ts, T, netGraph, dkfAlpha, dkfBeta, dkfDelta, b, P0);
  ckf     = CKF(plant, Ts, T);

  bLocal = rdkfloc.b;
  bSens  = bLocal(netGraph.Nodes.isSensor);
  fprintf('Local tolerances b^i: min=%.3g median=%.3g max=%.3g (global b=%.3g)\n', ...
          min(bSens), median(bSens), max(bSens), b);

  ckfRmse     = zeros(totalRuns, T + 1);
  crkfRmse    = zeros(totalRuns, T + 1);
  dkfRmse     = zeros(totalRuns, T + 1);
  rdkfRmse    = zeros(totalRuns, T + 1);
  rdkflocRmse = zeros(totalRuns, T + 1);
  dkfTxRate     = zeros(totalRuns, T + 1);
  rdkfTxRate    = zeros(totalRuns, T + 1);
  rdkflocTxRate = zeros(totalRuns, T + 1);

  tic
  parfor r = 1:totalRuns
    % Deterministic per-run trajectory => identical data for every b and filter.
    rng(trajSeedBase + r);
    s = lfm.simulate(x0);
    X = s.X; Y = s.Y;

    f = ckf.estimate(x0_hat, P0, X, Y);      ckfRmse(r, :) = f.RMSE;
    f = crkf.estimate(x0_hat, P0, X, Y);     crkfRmse(r, :) = f.RMSE;
    f = dkf.estimate(x0_hat, P0, X, Y);
    dkfRmse(r, :) = f.RMSE;     dkfTxRate(r, :) = f.txRate;
    f = rdkf.estimate(x0_hat, P0, X, Y);
    rdkfRmse(r, :) = f.RMSE;    rdkfTxRate(r, :) = f.txRate;
    f = rdkfloc.estimate(x0_hat, P0, X, Y);
    rdkflocRmse(r, :) = f.RMSE; rdkflocTxRate(r, :) = f.txRate;
  end
  fprintf('%d runs in %.1f s\n', totalRuns, toc);

  % --- paired summary (the statistic that matters at these margins) -----------
  ssCkf  = mean(ckfRmse(:, w), 2);
  ssCrkf = mean(crkfRmse(:, w), 2);
  ssDkf  = mean(dkfRmse(:, w), 2);
  ssRdkf = mean(rdkfRmse(:, w), 2);
  ssLoc  = mean(rdkflocRmse(:, w), 2);

  dCen = ssCkf - ssCrkf;                 % >0: CRKF better than CKF
  dDis = ssDkf - ssRdkf;                 % >0: RDKF better than DKF
  seCen = std(dCen) / sqrt(totalRuns);
  seDis = std(dDis) / sqrt(totalRuns);

  summary(ib, :) = [b, mean(ssDkf), mean(ssRdkf), ...
                    100 * mean(dDis) / mean(ssDkf), mean(dDis) / seDis, ...
                    100 * mean(dCen) / mean(ssCkf), mean(dCen) / seCen, ...
                    mean(mean(rdkfTxRate(:, w), 2)) / mean(mean(dkfTxRate(:, w), 2))];

  fprintf(['  CKF %.4f | CRKF %.4f (%+.2f%%, %.1f SE)\n' ...
           '  DKF %.4f | RDKF %.4f (%+.2f%%, %.1f SE) | RDKFLOC %.4f\n' ...
           '  TX: DKF %.4f  RDKF %.4f  (ratio %.2f)\n'], ...
          mean(ssCkf), mean(ssCrkf), 100 * mean(dCen) / mean(ssCkf), mean(dCen) / seCen, ...
          mean(ssDkf), mean(ssRdkf), 100 * mean(dDis) / mean(ssDkf), mean(dDis) / seDis, ...
          mean(ssLoc), mean(mean(dkfTxRate(:, w), 2)), mean(mean(rdkfTxRate(:, w), 2)), ...
          mean(mean(rdkfTxRate(:, w), 2)) / mean(mean(dkfTxRate(:, w), 2)));

  %% Save this b
  results = struct( ...
    'ckfRmse', ckfRmse, 'crkfRmse', crkfRmse, 'dkfRmse', dkfRmse, ...
    'rdkfRmse', rdkfRmse, 'rdkflocRmse', rdkflocRmse, ...
    'dkfTxRate', dkfTxRate, 'rdkfTxRate', rdkfTxRate, 'rdkflocTxRate', rdkflocTxRate);
  extras = struct('totalRuns', totalRuns, 'bLocal', bLocal, ...
                  'bGrid', bGrid, 'trajSeedBase', trajSeedBase, ...
                  'pairedDistributed', dDis, 'pairedCentralized', dCen, ...
                  'ssWindow', [w(1) w(end)]);

  % collectParams re-sources sst3dParams, so the overrides must be re-applied or
  % the saved file (and its name) would describe a run that did not happen.
  params = collectParams('sst3dParams');
  params.sensorCount   = sensorCount;
  params.dkfDelta      = dkfDelta;
  params.klTolerance   = b;
  params.lfmKlTolerance = lfmKlTolerance;
  params.T             = T;
  params.totalRuns     = totalRuns;
  params.bGrid         = bGrid;

  savedPaths(ib) = string(saveRun(mfilename, params, extras, netGraph, results, struct()));
end

%% Sweep summary
fprintf('\n================ SWEEP SUMMARY (%d paired runs) ================\n', totalRuns);
fprintf('%8s %9s %9s %10s %8s %11s %8s %8s %s\n', 'b', 'DKF', 'RDKF', ...
        'RDKF gain', 'SE', 'CRKF gain', 'SE', 'TXratio', 'trigger');
for ib = 1:numel(bGrid)
  s = summary(ib, :);
  fprintf('%8g %9.4f %9.4f %9.2f%% %8.1f %10.2f%% %8.1f %8.2f  %s\n', ...
          s(1), s(2), s(3), s(4), s(5), s(6), s(7), s(8), ...
          ternaryStr(s(1) <= 0.0088, 'feasible (F)', 'INFEASIBLE (F)'));
end
fprintf(['\n"RDKF gain" is the paired improvement over DKF; "CRKF gain" the paired\n' ...
         'improvement of the centralized robust filter over CKF at the SAME b.\n' ...
         'Condition (F) predicts the RDKF gain turns negative above b = 0.0088.\n']);
for ib = 1:numel(bGrid)
  fprintf('  %s\n', savedPaths(ib));
end

function s = ternaryStr(c, a, b)
  if c, s = a; else, s = b; end
end
