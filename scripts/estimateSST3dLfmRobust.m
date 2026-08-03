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
ssRdkfByB = zeros(totalRuns, numel(bGrid));   % kept for the iso-TX arm below
ssDkfRef  = zeros(totalRuns, 1);
txRdkfByB = zeros(numel(bGrid), 1);
txDkfRef  = 0;

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

  % DKF does not depend on b, so its column is identical every pass -- keeping it
  % once is enough, and overwriting it each pass is a free consistency check.
  ssRdkfByB(:, ib) = ssRdkf;
  txRdkfByB(ib)    = mean(mean(rdkfTxRate(:, w), 2));
  ssDkfRef         = ssDkf;
  txDkfRef         = mean(mean(dkfTxRate(:, w), 2));
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
  params = collectParams();
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

%% Iso-TX arm: is RDKF's gain bought with bandwidth?
%
% The b sweep above shows RDKF beating DKF, but at a higher transmission rate
% (ratio 1.25-1.31 at feasible b). That is on-model -- the paper's own Fig. 3
% reports 0.38 vs 0.28 = 1.36x -- but it is not a clean comparison: some of RDKF's
% RMSE advantage may simply be the extra packets, not the robustness.
%
% The test is to spend DKF's OWN trigger budget up to RDKF's transmission rate and
% re-run the comparison. In eq. (9) silence needs BOTH ||e||^2_Omega <= alpha and
% the Loewner sandwich, so DKF transmits more when alpha falls (innovation gate
% tightens) or when beta falls (lower Loewner bound tightens). Sweeping both spans
% DKF's whole (TX, RMSE) trade-off curve.
%
% Two readings, and the frontier is the stronger one:
%
%   * matched point -- the DKF config whose measured TX is closest to RDKF's. A
%     direct paired head-to-head at equal bandwidth.
%   * full frontier -- if RDKF's (TX, ssRMSE) point lies BELOW DKF's Pareto
%     frontier, the gain is real robustness. If it lies on or above it, RDKF is
%     just buying accuracy with bandwidth and the honest conclusion is that the
%     paper's Fig. 3 comparison is confounded.
%
% Held at b = isoB, the best feasible row of the sweep. Same network, plant, LFM
% and trajectory seeds, so every number here is paired with the sweep above.

% Grid sized from a T=250 probe of the same geometry (scratchpad lever.m), which
% found alpha to be a WEAK lever and beta a strong one: at alpha=10 the TX runs
% 0.415 / 0.522 / 0.706 / 1.000 for beta = 0.2 / 0.1 / 0.05 / 0.02. The grid below
% spans TX ~0.41 to ~0.97, so RDKF's 0.54 is bracketed from both sides rather than
% approached from below. The probe also showed DKF's ssRMSE *improving* with
% bandwidth (1.373 -> 1.327 at TX 0.528), which is exactly why this arm is needed.
isoB = 0.005;
isoAlphaGrid = [10 1 0.3 0.1 0.01];   % lower alpha => tighter innovation gate
isoBetaGrid  = [0.2 0.1 0.05];        % lower beta  => tighter Loewner bound

ibIso = find(bGrid == isoB, 1);
assert(~isempty(ibIso), 'isoB = %g must be one of bGrid', isoB);
ssRdkfIso = ssRdkfByB(:, ibIso);
txRdkfIso = txRdkfByB(ibIso);

fprintf(['\n\n================ ISO-TX ARM (b = %g) ================\n' ...
         'Target: DKF at RDKF''s transmission rate %.4f (RDKF ssRMSE %.4f).\n' ...
         'Baseline DKF (alpha=%g, beta=%g): TX %.4f, ssRMSE %.4f\n' ...
         'Sweeping alpha %s x beta %s over %d paired runs each.\n'], ...
        isoB, txRdkfIso, mean(ssRdkfIso), dkfAlpha, dkfBeta, txDkfRef, ...
        mean(ssDkfRef), mat2str(isoAlphaGrid), mat2str(isoBetaGrid), totalRuns);

nA = numel(isoAlphaGrid);
nB2 = numel(isoBetaGrid);
isoSs = zeros(totalRuns, nA, nB2);
isoTx = zeros(totalRuns, nA, nB2);

tic
parfor r = 1:totalRuns
  rng(trajSeedBase + r);            % same trajectory as the corresponding sweep run
  s = lfm.simulate(x0);
  X = s.X; Y = s.Y;

  ss = zeros(nA, nB2);
  tx = zeros(nA, nB2);
  for ia = 1:nA
    for ibb = 1:nB2
      f = DKF(plant, Ts, T, netGraph, isoAlphaGrid(ia), isoBetaGrid(ibb), dkfDelta) ...
            .estimate(x0_hat, P0, X, Y);
      ss(ia, ibb) = mean(f.RMSE(w));
      tx(ia, ibb) = mean(f.txRate(w));
    end
  end
  isoSs(r, :, :) = ss;
  isoTx(r, :, :) = tx;
end
fprintf('%d configs x %d runs in %.1f s\n', nA * nB2, totalRuns, toc);

%% Frontier table and the matched-TX verdict
fprintf('\n%8s %8s %9s %10s %12s %10s %s\n', 'alpha', 'beta', 'DKF TX', ...
        'DKF ssRMSE', 'vs RDKF', 'paired SE', 'who wins');
bestGap = inf; bestIdx = [0 0];
for ia = 1:nA
  for ibb = 1:nB2
    tx = mean(isoTx(:, ia, ibb));
    d = isoSs(:, ia, ibb) - ssRdkfIso;        % >0 means RDKF is still better
    se = std(d) / sqrt(totalRuns);
    rel = 100 * mean(d) / mean(isoSs(:, ia, ibb));
    tstat = mean(d) / (se + eps);
    if tstat > 2,       who = 'RDKF';
    elseif tstat < -2,  who = 'DKF';
    else,               who = 'tie';
    end
    fprintf('%8g %8g %9.4f %10.4f %11.2f%% %10.1f  %s\n', ...
            isoAlphaGrid(ia), isoBetaGrid(ibb), tx, mean(isoSs(:, ia, ibb)), ...
            rel, tstat, who);
    if abs(tx - txRdkfIso) < bestGap
      bestGap = abs(tx - txRdkfIso); bestIdx = [ia ibb];
    end
  end
end

ia = bestIdx(1); ibb = bestIdx(2);
dMatch = isoSs(:, ia, ibb) - ssRdkfIso;
seMatch = std(dMatch) / sqrt(totalRuns);
fprintf(['\n--- MATCHED-TX HEAD-TO-HEAD ---\n' ...
         'Closest DKF config: alpha=%g, beta=%g -> TX %.4f (RDKF TX %.4f, gap %.4f)\n' ...
         '  DKF  ssRMSE %.4f\n  RDKF ssRMSE %.4f  (b = %g)\n' ...
         '  paired difference %+.4f (%+.2f%%), %.1f SE\n'], ...
        isoAlphaGrid(ia), isoBetaGrid(ibb), mean(isoTx(:, ia, ibb)), txRdkfIso, ...
        bestGap, mean(isoSs(:, ia, ibb)), mean(ssRdkfIso), isoB, ...
        mean(dMatch), 100 * mean(dMatch) / mean(isoSs(:, ia, ibb)), ...
        mean(dMatch) / (seMatch + eps));
% The stricter test: the BEST DKF config that spends no more bandwidth than RDKF.
% Matching TX exactly is generous to RDKF, since the nearest config may sit above
% its TX. This asks whether RDKF beats every DKF tuning that is at least as frugal.
txMeans = squeeze(mean(isoTx, 1));
ssMeans = squeeze(mean(isoSs, 1));
frugal = txMeans <= txRdkfIso + 1e-9;
if ~any(frugal(:))
  fprintf(['\n--- FRUGAL-DKF TEST: no DKF config in the grid transmits as little as\n' ...
           'RDKF (%.4f). Widen isoAlphaGrid/isoBetaGrid upward.\n'], txRdkfIso);
else
  cand = ssMeans; cand(~frugal) = inf;
  [~, k] = min(cand(:));
  [ja, jb] = ind2sub(size(cand), k);
  dFrug = isoSs(:, ja, jb) - ssRdkfIso;
  seFrug = std(dFrug) / sqrt(totalRuns);
  fprintf(['\n--- FRUGAL-DKF TEST (best DKF with TX <= RDKF''s %.4f) ---\n' ...
           'alpha=%g, beta=%g -> TX %.4f, ssRMSE %.4f\n' ...
           'RDKF ssRMSE %.4f at TX %.4f\n' ...
           '  paired difference %+.4f (%+.2f%%), %.1f SE  => %s\n'], ...
          txRdkfIso, isoAlphaGrid(ja), isoBetaGrid(jb), txMeans(ja, jb), ...
          ssMeans(ja, jb), mean(ssRdkfIso), txRdkfIso, mean(dFrug), ...
          100 * mean(dFrug) / ssMeans(ja, jb), mean(dFrug) / (seFrug + eps), ...
          ternaryStr(mean(dFrug) / (seFrug + eps) > 2, 'RDKF survives', ...
                     'RDKF gain is bandwidth-confounded'));
end

fprintf(['\nRead this against the frontier: a positive difference at |SE| > 2 in\n' ...
         'BOTH tests means RDKF''s advantage survives equal bandwidth and is genuine\n' ...
         'robustness. A tie or negative in the frugal test means the sweep''s gain\n' ...
         'was bought with packets, and the paper''s Fig. 3 comparison (RDKF 0.38 vs\n' ...
         'DKF 0.28) is bandwidth-confounded -- a finding in its own right.\n']);

isoExtras = struct('isoB', isoB, 'isoAlphaGrid', isoAlphaGrid, ...
                   'isoBetaGrid', isoBetaGrid, 'isoSs', isoSs, 'isoTx', isoTx, ...
                   'ssRdkfIso', ssRdkfIso, 'txRdkfIso', txRdkfIso, ...
                   'ssDkfRef', ssDkfRef, 'txDkfRef', txDkfRef, ...
                   'totalRuns', totalRuns, 'trajSeedBase', trajSeedBase, ...
                   'ssWindow', [w(1) w(end)]);
isoParams = collectParams();
isoParams.sensorCount    = sensorCount;
isoParams.dkfDelta       = dkfDelta;
isoParams.klTolerance    = isoB;
isoParams.lfmKlTolerance = lfmKlTolerance;
isoParams.T              = T;
isoParams.totalRuns      = totalRuns;
fprintf('\n%s\n', saveRun([mfilename '_isoTx'], isoParams, isoExtras, netGraph, ...
                          struct('isoSs', isoSs, 'isoTx', isoTx), struct()));

function s = ternaryStr(c, a, b)
  if c, s = a; else, s = b; end
end
