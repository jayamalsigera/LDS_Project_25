%% Fusion-weight unit test — Pi must be row-stochastic over the set fusion() visits
%
% Regression test for the edge-direction bug described in
% docs/rdkf-tx-saturation-analysis.md §6.4.1.
%
% Every distributed filter in this repo fuses over `inedges(G, i)`, i.e. over the
% in-neighborhood N_i of node i:
%
%   for i = 1:N
%     [~, nids] = inedges(self.G, i);
%     for j = nids'
%       ... + self.Pi(i, j) * (information from j)
%
% For that fusion to be a convex combination, `calcFusionWeights` must build Pi
% over the SAME set, so that
%
%   sum_{j in {i} union N_i} Pi(i, j) == 1     for every i,
%
% and Pi must place no mass on nodes outside that set (any such mass is silently
% dropped by the fusion loop, which never visits them).
%
% `calcFusionWeights` originally built Pi over the OUT-neighborhood. On a
% symmetric graph the two coincide, so every graph the repo builds
% (`createSpatialNetwork` makes an undirected `graph` and converts it) hid the
% bug. On a genuine digraph the row sums over the visited set fell to a median of
% ~0.2, so ~80% of the network information leaked away each iteration, Omega
% decayed geometrically to zero, and `findTheta` ended up inverting a singular
% matrix — with no error raised anywhere.
%
% This test therefore checks BOTH a symmetric graph (must still pass) and a
% deliberately asymmetric digraph (the case that used to fail).

clear; clc;
fprintf('=== FUSION-WEIGHT UNIT TEST ===\n\n');

TOL = 1e-12;
rng(42);

cases = {};

% --- Case 1: small hand-built digraph with deliberately asymmetric edges ------
% Node 1 -> 2, 2 -> 3, 3 -> 1 (a cycle: every in-degree is 1), plus the extra
% one-way edges 1 -> 3 and 1 -> 4, and 4 -> 1. Self-loops on all nodes.
% In-degrees (excl. self): node1: {3,4} = 2, node2: {1} = 1, node3: {2,1} = 2,
% node4: {1} = 1. Out-degrees differ (node1 has 3), so out- vs in-neighborhood
% weighting give different answers here.
Ah = false(4, 4);
Ah(1, 2) = true; Ah(2, 3) = true; Ah(3, 1) = true;
Ah(1, 3) = true; Ah(1, 4) = true; Ah(4, 1) = true;
Ah(1:5:end) = true;                                    % self-loops
cases{end+1} = struct('name', 'hand-built asymmetric digraph (N=4)', ...
                      'G', digraph(sparse(Ah)));

% --- Case 2: random Erdos-Renyi digraph at the paper's 4% connection density --
% This is the topology of Ghion & Zorzi (2023) Section 6. It is genuinely
% directed: rand(N,N) < p is asymmetric with probability ~1.
N = 60;
Aer = rand(N, N) < 0.04;
Aer(sub2ind([N N], 1:N, [2:N 1])) = true;              % Hamiltonian cycle: strong connectivity
Aer(1:N+1:end) = true;                                 % self-loops
cases{end+1} = struct('name', sprintf('ER 4%% digraph (N=%d)', N), ...
                      'G', digraph(sparse(Aer)));

% --- Case 3: the shipped spatial network (symmetric) — must not regress -------
cases{end+1} = struct('name', 'createSpatialNetwork (symmetric, N=40 S=30)', ...
                      'G', createSpatialNetwork(40, 30, 5000));

% Both weight builders must satisfy the same contract: RDKF/SRDKF use
% calcFusionWeights, DKF/SDKF/DSEACP use calcMetropolisWeights, and every one of
% them fuses over inedges.
builders = {@calcFusionWeights, @calcMetropolisWeights};
bnames   = {'calcFusionWeights', 'calcMetropolisWeights'};

allPass = true;

for bi = 1:numel(builders)
fprintf('---- %s ----\n', bnames{bi});
for c = 1:numel(cases)
  G = cases{c}.G;
  Nn = numnodes(G);
  Pi = builders{bi}(G);

  % Row sums over exactly the set the estimators' fusion loop iterates.
  visitedSum = zeros(Nn, 1);
  strayMass  = zeros(Nn, 1);
  for i = 1:Nn
    [~, nids] = inedges(G, i);
    visited = unique([i; nids(:)]);
    visitedSum(i) = sum(Pi(i, visited));

    mask = true(1, Nn);
    mask(visited) = false;
    strayMass(i) = sum(abs(Pi(i, mask)));
  end

  % A non-negative, convex combination is also required for the fused Omega to
  % stay positive definite.
  nonNeg = all(Pi(:) >= -TOL);

  okSum   = all(abs(visitedSum - 1) < TOL);
  okStray = all(strayMass < TOL);
  pass    = okSum && okStray && nonNeg;
  allPass = allPass && pass;

  isSym = issymmetric(double(adjacency(G) > 0));
  fprintf('%-48s %s\n', cases{c}.name, tf(pass));
  fprintf('  symmetric adjacency          : %d\n', isSym);
  fprintf('  row sum over visited set     : min %.6f  median %.6f  max %.6f   %s\n', ...
          min(visitedSum), median(visitedSum), max(visitedSum), tf(okSum));
  fprintf('  mass on never-visited nodes  : max %.3e   %s\n', max(strayMass), tf(okStray));
  fprintf('  all weights non-negative     : %s\n\n', tf(nonNeg));
end
end

% Cross-check: the fused information of a node whose in-neighbors all carry the
% same Omega must equal that Omega exactly (convexity), and it must NOT shrink.
% This is the property whose violation caused the geometric decay to zero.
G  = cases{2}.G;
Nn = numnodes(G);
Pi = calcFusionWeights(G);
Om = diag([2.6 4.3 5.5 46.6 53.1 59.0]);              % steady-state-like Omega
worst = 0;
for i = 1:Nn
  [~, nids] = inedges(G, i);
  fused = zeros(6);
  for j = unique([i; nids(:)])'
    fused = fused + Pi(i, j) * Om;
  end
  worst = max(worst, norm(fused - Om) / norm(Om));
end
okFuse  = worst < TOL;
allPass = allPass && okFuse;
fprintf('identical-neighbor fusion is idempotent: worst rel. error %.3e   %s\n\n', ...
        worst, tf(okFuse));

if allPass
  fprintf('=== FUSION-WEIGHT UNIT TEST PASSED ===\n');
else
  error('fusionWeightsUnitTest: FAILED — see the per-case report above.');
end

function s = tf(b)
  if b, s = 'PASS'; else, s = 'FAIL'; end
end
