%% Validation report for the per-node local tolerances b^i (RDKFLOC/SRDKFLOC)
%
% Standalone printout of the local tolerances computed by
% computeLocalTolerances from the global least-favorable model. No estimation,
% no Monte Carlo -- just the offline b^i and the sanity checks that justify
% Algorithm 2 (Ghion & Zorzi 2023, Section 5 / Fig. 4).
%
% Headline check: the median b^i should sit near the ~1e-3 empirical tuning
% optimum, and every sensor b^i must satisfy 0 < b^i <= global tolerance.

clear;
clc;
rng(42);

%% Parameters
sst2dParams;

%% Network and plant
netGraph = createSpatialNetwork(nodeCount, sensorCount, maxLength);
plant    = SingleTarget2dModel(Ts, sensorCount, outputNoiseStd, T, turnRate);

%% Compute local tolerances
fprintf('Computing local tolerances from global LFM (b = %.3g)...\n', lfmKlTolerance);
bLocal = computeLocalTolerances(plant, P0, lfmKlTolerance, netGraph);

isSensor = netGraph.Nodes.isSensor;
bSens    = bLocal(isSensor);

%% Summary
fprintf('\n=== Local tolerances b^i ===\n');
fprintf('Sensor nodes            : %d\n', numel(bSens));
fprintf('min / median / max b^i  : %.4g / %.4g / %.4g\n', ...
        min(bSens), median(bSens), max(bSens));
fprintf('global tolerance b      : %.4g\n', lfmKlTolerance);
fprintf('all 0 < b^i <= b        : %s\n', ...
        mat2str(all(bSens > 0 & bSens <= lfmKlTolerance * (1 + 1e-6))));

% Per-axis groups: sensor nodes 1..S/2 measure p_x, S/2+1..S measure p_y. The
% two medians coincide because the plant is symmetric under x<->y (identical
% dynamics and noise per axis), so the homogeneous network yields a single b^i.
S       = sensorCount;
axisPx  = bLocal(1:floor(S / 2));
axisPy  = bLocal(floor(S / 2) + 1:S);
fprintf('\n=== Per-axis clusters ===\n');
fprintf('p_x sensors  median b^i : %.4g  (n=%d)\n', median(axisPx), numel(axisPx));
fprintf('p_y sensors  median b^i : %.4g  (n=%d)\n', median(axisPy), numel(axisPy));

fprintf('\nHeadline check: median b^i (%.4g) vs empirical tuning optimum (~1e-3).\n', ...
        median(bSens));
