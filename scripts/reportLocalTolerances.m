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
sst3dParams;

%% Network and plant
netGraph = createSpatialNetwork(nodeCount, sensorCount, maxLength);
plant    = SingleTarget3dModel(Ts, sensorCount, noiseScale, T);

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

% Per-type groups: SingleTarget3dModel splits the sensors into three contiguous
% blocks by measured coordinate pair (px/py, py/pz, px/pz), sized as evenly as S
% allows. The medians coincide because the plant is symmetric under coordinate
% permutation (identical dynamics and noise per axis), so the homogeneous network
% yields a single b^i.
S        = sensorCount;
nPerType = floor(S / 3) * ones(1, 3);
nPerType(1:(S - sum(nPerType))) = nPerType(1:(S - sum(nPerType))) + 1;
typeEnd   = cumsum(nPerType);
typeStart = [1, typeEnd(1:2) + 1];
typeName  = {'px,py', 'py,pz', 'px,pz'};
fprintf('\n=== Per-type clusters ===\n');
for t = 1:3
  bType = bLocal(typeStart(t):typeEnd(t));
  fprintf('%s sensors  median b^i : %.4g  (n=%d)\n', ...
          typeName{t}, median(bType), numel(bType));
end

fprintf('\nHeadline check: median b^i (%.4g) vs empirical tuning optimum (~1e-3).\n', ...
        median(bSens));
